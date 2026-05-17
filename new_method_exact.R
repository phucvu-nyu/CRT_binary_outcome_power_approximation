suppressPackageStartupMessages({
  # all used functions are in base/stats, so no extra packages required
})

############################################################
# 0. Generate betas from an effect map
############################################################
# effect_map is a named list, e.g.:
# list("component_1" = 0.25,
#      "component_1*component_2" = 0.125, ...)
#
# Betas are ordered as:
#  - main effects: component_1, ..., component_k
#  - interactions: all 2-way, then 3-way, ..., up to interaction_order
#    using combn() order, matching compute_factorial_probabilities().
generate_betas <- function(k_factor = 5, interaction_order = 5, effect_map = NULL) {
  num_main <- k_factor
  num_interactions <- if (interaction_order >= 2) {
    sum(choose(k_factor, 2:interaction_order))
  } else 0
  total_betas <- num_main + num_interactions
  
  betas <- numeric(total_betas)
  component_names <- paste0("component_", 1:k_factor)
  
  # Build the interaction name list in the same order as the design matrix
  interaction_names <- character(0)
  if (interaction_order >= 2) {
    for (order in 2:interaction_order) {
      combos <- apply(combn(component_names, order), 2, function(x) paste(x, collapse = "*"))
      interaction_names <- c(interaction_names, combos)
    }
  }
  
  # Lookup for interactions: offset by k_factor for main effects
  if (length(interaction_names) > 0) {
    lookup <- setNames(seq_along(interaction_names) + k_factor, interaction_names)
  } else {
    lookup <- named(integer(0))
  }
  
  if (!is.null(effect_map)) {
    for (term in names(effect_map)) {
      if (term %in% component_names) {
        idx <- which(component_names == term)
      } else if (term %in% names(lookup)) {
        idx <- lookup[term]
      } else {
        stop(paste("Term not found in main effects or interactions:", term))
      }
      betas[idx] <- effect_map[[term]]
      cat(term, "-> beta =", effect_map[[term]], "at index", idx, "\n")
    }
  }
  
  betas
}

############################################################
# 1. Factorial truth generator (fixed-effect structure)
############################################################
# Given betas and baseline probability at the all -1 corner,
# build:
#  - design_matrix: main-effects ±1 for each component
#  - design_full: main + interactions up to interaction_order
#  - eta_all: linear predictors (no random effects)
compute_factorial_probabilities <- function(
    p_baseline,
    betas,
    k_factor = 5,
    interaction_order = 5,
    design_matrix = NULL
) {
  # 1) Full factorial main-effects design if not provided
  if (is.null(design_matrix)) {
    design_matrix <- as.matrix(expand.grid(rep(list(c(-1, 1)), k_factor)))
    colnames(design_matrix) <- paste0("component_", 1:k_factor)
  } else {
    design_matrix <- as.matrix(design_matrix)
  }
  
  # 2) Build design matrix with interactions up to interaction_order
  X <- design_matrix
  design_full <- X
  term_names  <- colnames(X)  # main effects first
  
  if (interaction_order >= 2) {
    for (order in 2:interaction_order) {
      combos <- combn(colnames(X), order, simplify = FALSE)
      for (cmb in combos) {
        term <- apply(X[, cmb, drop = FALSE], 1, prod)
        design_full <- cbind(design_full, term)
        term_names  <- c(term_names, paste(cmb, collapse = "*"))
      }
    }
  }
  
  # 3) Calibrate intercept so that all components = -1 gives p_baseline
  x_minus1 <- rep(-1, k_factor)
  interaction_minus1 <- numeric(0)
  
  if (interaction_order >= 2) {
    for (order in 2:interaction_order) {
      combos <- combn(1:k_factor, order)
      for (j in seq_len(ncol(combos))) {
        interaction_minus1 <- c(interaction_minus1, prod(x_minus1[combos[, j]]))
      }
    }
  }
  
  signs_full <- c(x_minus1, interaction_minus1)
  if (length(signs_full) != length(betas)) {
    stop("Length of betas does not match number of terms implied by k_factor and interaction_order.")
  }
  
  eta_minus1 <- sum(betas * signs_full)
  beta0 <- qlogis(p_baseline) - eta_minus1
  
  # 4) Linear predictors for all design cells
  eta_all <- as.numeric(beta0 + design_full %*% betas)
  p_all   <- plogis(eta_all)
  
  list(
    overall_avg       = mean(p_all),
    design_matrix     = design_matrix,   # main-effects ±1
    design_full       = design_full,     # main + interactions
    term_names        = term_names,
    eta_all           = eta_all,
    p_all_conditional = p_all
  )
}

############################################################
# 2. Balanced sampling of design rows for K clusters
############################################################
balanced_sample_indices <- function(J, K) {
  base <- rep(K %/% J, J)
  remainder <- K - sum(base)
  if (remainder > 0) {
    extra <- sample(seq_len(J), remainder, replace = FALSE)
    base[extra] <- base[extra] + 1L
  }
  sample(rep(seq_len(J), base), size = K, replace = FALSE)
}

expand_design_for_clusters <- function(X_full, w_full, K) {
  J <- nrow(X_full)
  idx <- balanced_sample_indices(J, K)
  X_cluster <- X_full[idx, , drop = FALSE]
  w_cluster <- w_full[idx]
  list(X_cluster = X_cluster, w_cluster = w_cluster)
}

compute_covariance_from_w <- function(X_cluster, w_cluster) {
  W_inv <- diag(1 / w_cluster)
  solve(t(X_cluster) %*% W_inv %*% X_cluster)
}

############################################################
# 3. General cluster size generator
############################################################
# dist = "uniform":
#   - if min_m/max_m not provided, use [round(m/2), round(3m/2)]
# dist = "gaussian":
#   - need sd_m
#   - sample N(mean_m, sd_m), round, floor at 1, optionally cap at max_m
sample_cluster_sizes_general <- function(
    mean_m,
    B,
    dist = c("uniform", "gaussian"),
    min_m = NULL,
    max_m = NULL,
    sd_m  = NULL
) {
  dist <- match.arg(dist)
  
  if (dist == "uniform") {
    if (is.null(min_m)) min_m <- round(mean_m / 2)
    if (is.null(max_m)) max_m <- round(3 * mean_m / 2)
    
    min_m <- max(1L, as.integer(min_m))
    max_m <- max(min_m, as.integer(max_m))
    
    return(sample(seq(min_m, max_m), size = B, replace = TRUE))
  }
  
  if (dist == "gaussian") {
    if (is.null(sd_m)) stop("For dist = 'gaussian', sd_m must be provided.")
    
    sizes <- round(rnorm(B, mean = mean_m, sd = sd_m))
    sizes[sizes < 1] <- 1L
    
    if (!is.null(max_m)) {
      sizes[sizes > max_m] <- max_m
    }
    return(as.integer(sizes))
  }
}

############################################################
# 4. Var[logit(cluster mean)] via GLMM Monte Carlo (one cell)
############################################################
# For a single factorial cell j with linear predictor eta_j and
# random intercept variance tau2:
#  - draw cluster sizes m ~ chosen distribution
#  - draw b ~ N(0, tau2)
#  - draw Y_ij ~ Bernoulli(logit^{-1}(eta_j + b))
#  - compute cluster means and their logits
#  - return empirical Var[logit(cluster mean)]
var_logit_cluster_mc_glmm <- function(
    eta_j,
    tau2,
    mean_m,
    dist = "uniform",
    min_m = NULL,
    max_m = NULL,
    sd_m  = NULL,
    B_mc  = 5000,
    eps   = 0.5
) {
  # cluster sizes
  m_vec <- sample_cluster_sizes_general(
    mean_m = mean_m,
    B      = B_mc,
    dist   = dist,
    min_m  = min_m,
    max_m  = max_m,
    sd_m   = sd_m
  )
  
  # random intercepts
  b_vec <- rnorm(B_mc, mean = 0, sd = sqrt(tau2))
  
  # cluster-specific probabilities
  p_vec <- plogis(eta_j + b_vec)
  
  # binomial totals and proportions
  k_vec <- rbinom(B_mc, size = m_vec, prob = p_vec)
  # continuity correction to avoid logit(0) / logit(1)
  ybar_vec <- (k_vec + eps) / (m_vec + 2 * eps)
  
  var(qlogis(ybar_vec))
}

############################################################
# 5. W = diag(Var[logit(cluster mean)]) via GLMM MC
############################################################
compute_W_glmm_mc <- function(
    eta_all,
    tau2,
    mean_m,
    dist = "uniform",
    min_m = NULL,
    max_m = NULL,
    sd_m  = NULL,
    B_mc  = 5000,
    eps   = 0.5
) {
  J <- length(eta_all)
  w <- numeric(J)
  
  for (j in seq_len(J)) {
    w[j] <- var_logit_cluster_mc_glmm(
      eta_j = eta_all[j],
      tau2  = tau2,
      mean_m = mean_m,
      dist   = dist,
      min_m  = min_m,
      max_m  = max_m,
      sd_m   = sd_m,
      B_mc   = B_mc,
      eps    = eps
    )
  }
  
  list(
    W = diag(w),
    w = w
  )
}

############################################################
# 6. Main function: GLMM-based MC analytic power
############################################################
# This:
#  - builds the full factorial truth (fixed effects)
#  - computes tau2 from ICC on latent scale
#  - uses GLMM MC to get Var[logit(cluster mean)] for each cell
#  - uses GLS approximation + noncentral t to get power for component_1
estimate_power_replicated_glmm_mc <- function(
    K = 50,               # number of clusters
    B = 100,              # replicates over random cluster allocation
    p_baseline = 0.40,    # baseline probability when all components = -1
    betas,
    k_factor = 5,
    fit_interaction_order = 2,  # order used in the fitted model
    icc_latent = 0.05,    # ICC on latent (logit) scale
    mean_m = 20,          # average cluster size
    dist  = "uniform",    # "uniform" or "gaussian"
    min_m = NULL,         # optional for uniform
    max_m = NULL,         # optional for uniform / gaussian cap
    sd_m  = NULL,         # required if dist = "gaussian"
    alpha = 0.05,
    seed  = 123,
    B_mc  = 5000          # MC reps per factorial cell for Var(logit(cluster mean))
) {
  set.seed(seed)
  
  ## STEP 1: Fixed-effects factorial truth
  test <- compute_factorial_probabilities(
    p_baseline        = p_baseline,
    betas             = betas,
    k_factor          = k_factor,
    interaction_order = k_factor
  )
  # test$eta_all: linear predictors
  # test$design_full: full design matrix (main + interactions)
  
  ## STEP 2: Random intercept variance from latent ICC
  logistic_var <- (pi^2) / 3
  tau2 <- logistic_var * icc_latent / (1 - icc_latent)
  
  ## STEP 3: W on logit scale via GLMM MC
  W_res <- compute_W_glmm_mc(
    eta_all = test$eta_all,
    tau2    = tau2,
    mean_m  = mean_m,
    dist    = dist,
    min_m   = min_m,
    max_m   = max_m,
    sd_m    = sd_m,
    B_mc    = B_mc
  )
  w_full <- W_res$w
  
  ## STEP 4: Full X for GLS (matches factorial terms)
  X_full <- cbind(1, test$design_full)
  colnames(X_full) <- c("(Intercept)", test$term_names)
  
  # contrast for component_1 (always column 2 in this coding)
  c_vec <- rep(0, ncol(X_full))
  c_vec[2] <- 1
  
  true_theta <- betas[1]  # main effect for component_1
  
  ## STEP 5: number of parameters in fitted model
  p_fit <- 1 + sum(choose(k_factor, 1:fit_interaction_order))
  
  ## STEP 6: simulate K-cluster allocations and compute power
  power_vec <- numeric(B)
  
  for (b in seq_len(B)) {
    cl <- expand_design_for_clusters(X_full, w_full, K)
    X_cluster <- cl$X_cluster
    w_cluster <- cl$w_cluster
    
    Sigma <- compute_covariance_from_w(X_cluster, w_cluster)
    
    df <- max(K - p_fit, 1)
    se_hat <- sqrt(as.numeric(t(c_vec) %*% Sigma %*% c_vec))
    lambda <- true_theta / se_hat
    tcrit  <- qt(1 - alpha/2, df)
    
    power_vec[b] <-
      pt(-tcrit, df = df, ncp = lambda) +
      (1 - pt(tcrit, df = df, ncp = lambda))
  }
  
  list(
    mean_power = mean(power_vec),
    sd_power   = sd(power_vec),
    all_power  = power_vec,
    df_used    = K - p_fit,
    p          = p_fit
  )
}

############################################################
# 7. Example usage (uncomment to run)
############################################################
# effect_map <- list(
#   "component_1" = 0.25,
#   "component_3" = 0,
#   "component_5" = 0,
#   "component_1*component_2" = 0,
#   "component_1*component_3" = 0,
#   "component_1*component_3*component_5" = 0,
#   "component_1*component_2*component_3*component_5" = 0
# )

# betas <- generate_betas(k_factor = 5, interaction_order = 5, effect_map = effect_map)
# 
# res <- estimate_power_replicated_glmm_mc(
#   K = 50,
#   B = 100,
#   p_baseline = 0.40,
#   betas = betas,
#   k_factor = 5,
#   fit_interaction_order = 2,
#   icc_latent = 0.1,
#   mean_m = 20,
#   dist = "uniform",   # or "gaussian"
#   alpha = 0.05,
#   B_mc = 5000
# )
# 
# print(res$mean_power)
