
# ---------------------------------------------------------
# Compute true average probabilities from factorial βs
# ---------------------------------------------------------
compute_factorial_probabilities <- function(
    p_baseline,
    betas,
    k_factor = 5,
    interaction_order = 5,
    design_matrix = NULL   # optional: pass fractional design rows
) {
  # 1. Construct the full factorial (+1/-1 for each component)
  if (is.null(design_matrix)) {
    # Full 2^k factorial
    design_matrix <- as.matrix(expand.grid(rep(list(c(-1, 1)), k_factor)))
    colnames(design_matrix) <- paste0("component_", 1:k_factor)
  } else {
    # Use provided fractional design
    design_matrix <- as.matrix(design_matrix)
  }
  
  # 2. Build interaction terms to match the betas vector
  X <- design_matrix
  design_full <- X
  term_names  <- colnames(X)   # main effects first
  
  if (interaction_order >= 2) {
    for (order in 2:interaction_order) {
      combos <- combn(colnames(X), order, simplify = FALSE)
      for (cmb in combos) {
        term <- apply(X[, cmb, drop = FALSE], 1, prod)
        term_name <- paste(cmb, collapse = "*")
        design_full <- cbind(design_full, term)
        term_names  <- c(term_names, term_name)
      }
    }
  }
  
  # 3. Insert betas in correct order and compute β0
  # Assume betas correspond exactly to the order produced by term_names
  # Calibrate intercept so that all components = -1 => p_baseline
  x_minus1 <- rep(-1, k_factor)
  interaction_minus1 <- c()
  
  if (interaction_order >= 2) {
    for (order in 2:interaction_order) {
      combos <- combn(1:k_factor, order)
      for (j in 1:ncol(combos)) {
        interaction_minus1 <- c(interaction_minus1, prod(x_minus1[combos[, j]]))
      }
    }
  }
  
  signs_full <- c(x_minus1, interaction_minus1)
  eta_minus1 <- sum(betas * signs_full)
  
  beta0 <- qlogis(p_baseline) - eta_minus1
  
  # 4. Compute linear predictors and conditional probabilities
  eta_all <- beta0 + as.numeric(design_full %*% betas)
  p_all  <- plogis(eta_all)   # conditional (b = 0), still useful to keep
  
  # 5. Average probabilities for each component = -1 or +1 (conditional)
  per_component <- list()
  for (j in 1:k_factor) {
    comp_name <- paste0("component_", j)
    per_component[[comp_name]] <- list(
      avg_minus1 = mean(p_all[design_matrix[, j] == -1]),
      avg_plus1  = mean(p_all[design_matrix[, j] == +1])
    )
  }
  
  results <- list(
    overall_avg     = mean(p_all),     # conditional avg
    per_component   = per_component,   # conditional per-component avgs
    design_matrix   = design_matrix,   # main-effects ±1 matrix
    design_full     = design_full,     # main + interactions (matches betas)
    term_names      = term_names,      # names for betas
    eta_all         = eta_all,         # fixed-effects linear predictors
    p_all_conditional = p_all          # conditional probs (b = 0)
  )
  
  return(results)
}


icc_marginal_from_pavg <- function(p_avg, icc_latent) {
  logistic_var <- (pi^2) / 3
  tau2 <- logistic_var * icc_latent / (1 - icc_latent)
  
  p_var <- p_avg * (1 - p_avg)
  
  V1 <- p_var
  V2 <- tau2 * p_var^2
  
  icc_marg <- V2 / (V1 + V2)
  
  return(list(
    tau2 = tau2,
    p_var = p_var,
    icc_marg = icc_marg
  ))
}

# ----------------------------------------------
# FUNCTION 1: Compute marginal mean under GLMM
# ----------------------------------------------
compute_marginal_mean <- function(eta, tau2, method = "GH", B = 1000) {
  eta  <- as.numeric(eta)
  J    <- length(eta)
  
  # Allow tau2 to be scalar or vector
  if (length(tau2) == 1L) {
    tau2 <- rep(tau2, J)
  } else if (length(tau2) != J) {
    stop("tau2 must be either length 1 or the same length as eta.")
  }
  
  # If all random-effect variances are (essentially) zero,
  # marginal mean = conditional mean
  if (all(tau2 < 1e-12)) {
    return(plogis(eta))
  }
  
  if (method == "MC") {
    # Monte Carlo integration, reusing the same standard normals
    b <- rnorm(B, mean = 0, sd = 1)
    # For each j, E[plogis(eta_j + sqrt(tau2_j) * b)]
    out <- vapply(
      seq_len(J),
      function(j) {
        mean(plogis(eta[j] + sqrt(tau2[j]) * b))
      },
      numeric(1)
    )
    return(out)
  }
  
  if (method == "GH") {
    # Gauss–Hermite quadrature with 20 nodes
    gh <- statmod::gauss.quad(20, kind = "hermite")
    base_nodes   <- gh$nodes
    weights      <- gh$weights / sqrt(pi)
    
    out <- vapply(
      seq_len(J),
      function(j) {
        nodes_j <- sqrt(2 * tau2[j]) * base_nodes
        sum(weights * plogis(eta[j] + nodes_j))
      },
      numeric(1)
    )
    return(out)
  }
  
  stop("Unknown method. Use 'GH' or 'MC'.")
}


compute_marginal_icc <- function(p, icc_latent) {
  logistic_var <- (pi^2) / 3
  tau2 <- logistic_var * icc_latent / (1 - icc_latent)
  
  p_var <- p * (1 - p)
  
  V1 <- p_var
  V2 <- tau2 * p_var^2
  
  icc_marg <- V2 / (V1 + V2)
  
  return(icc_marg)
}
# --------------------------------------------------
# FUNCTION 3 (CV version): Compute W using mean & SD
# --------------------------------------------------
compute_W <- function(p, rho, mean_m, sd_m) {
  p  <- as.numeric(p)
  rho <- as.numeric(rho)
  J  <- length(p)
  
  if (length(rho) == 1L) {
    rho <- rep(rho, J)
  } else if (length(rho) != J) {
    stop("rho must be length 1 or same length as p.")
  }
  
  if (mean_m <= 0) {
    stop("mean_m must be positive.")
  }
  
  # Coefficient of variation for cluster size
  cv <- if (sd_m < 1e-12) 0 else sd_m / mean_m
  
  m_bar <- mean_m
  
  # Unequal-size design effect adjustment:
  # replace (m_j - 1) by (m_bar * (1 + cv^2) - 1)
  m_star <- m_bar * (1 + cv^2)  # "effective" size in the DE term
  
  # 1. Variance on Y-scale for each cell (approximate)
  vY <- p * (1 - p) / m_bar * (1 + (m_star - 1) * rho)
  
  # 2. Delta-method variance on logit scale
  w <- vY / (p^2 * (1 - p)^2)
  
  # 3. Diagonal W
  W <- diag(w)
  
  return(list(
    W    = W,
    vY   = vY,
    w    = w,
    m_bar = m_bar,
    cv    = cv,
    m_star = m_star
  ))
}
compute_covariance <- function(X, W) {
  W_inv <- diag(1 / diag(W))
  Sigma <- solve(t(X) %*% W_inv %*% X)
  return(Sigma)
}
compute_power <- function(Sigma, c, true_theta, alpha = 0.05) {
  var_theta <- as.numeric(t(c) %*% Sigma %*% c)
  se_theta  <- sqrt(var_theta)
  
  lambda <- true_theta / se_theta
  zcrit  <- qnorm(1 - alpha/2)
  
  power <- pnorm(lambda - zcrit) + pnorm(-lambda - zcrit)
  return(power)
}
compute_power <- function(Sigma, c, true_theta, alpha = 0.05) {
  var_theta <- as.numeric(t(c) %*% Sigma %*% c)
  se_theta  <- sqrt(var_theta)
  
  lambda <- true_theta / se_theta
  zcrit  <- qnorm(1 - alpha/2)
  
  power <- pnorm(lambda - zcrit) + pnorm(-lambda - zcrit)
  return(power)
}

balanced_sample_indices <- function(J, K) {
  base <- rep(K %/% J, J)
  remainder <- K - sum(base)
  
  if (remainder > 0) {
    extra_arms <- sample(seq_len(J), remainder, replace = FALSE)
    base[extra_arms] <- base[extra_arms] + 1
  }
  
  alloc <- rep(seq_len(J), base)
  sample(alloc, size = K, replace = FALSE)
}
expand_design_for_clusters <- function(X_full, w_full, K) {
  J <- nrow(X_full)
  
  # Balanced sampling of cluster rows
  idx <- balanced_sample_indices(J, K)
  
  # Build cluster-level X and w
  X_cluster <- X_full[idx, , drop = FALSE]
  w_cluster <- w_full[idx]
  
  return(list(
    X_cluster = X_cluster,
    w_cluster = w_cluster
  ))
}
compute_covariance_from_w <- function(X_cluster, w_cluster) {
  W_inv <- diag(1 / w_cluster)
  Sigma <- solve(t(X_cluster) %*% W_inv %*% X_cluster)
  return(Sigma)
}
compute_power_small_sample <- function(Sigma, c, true_theta, alpha, K, p) {
  
  var_theta <- as.numeric(t(c) %*% Sigma %*% c)
  se_theta  <- sqrt(var_theta)
  
  lambda <- true_theta / se_theta
  
  # small-sample df
  df <- K - p   # p = number of parameters in model
  
  tcrit <- qt(1 - alpha/2, df)
  
  # two-sided power under noncentral t
  power <- pt(lambda - tcrit, df = df) + (1 - pt(lambda + tcrit, df = df))
  return(power)
}
compute_satterthwaite_df <- function(Sigma, c) {
  # numerator
  V <- as.numeric(t(c) %*% Sigma %*% c)
  
  # denominator = c^T (Sigma ∘ Sigma) c
  Sigma_sq <- Sigma * Sigma
  denom <- as.numeric(t(c) %*% Sigma_sq %*% c)
  
  df <- V^2 / denom
  return(df)
}

compute_power_satterthwaite <- function(Sigma, c, true_theta, alpha = 0.05) {
  # variance and se
  V <- as.numeric(t(c) %*% Sigma %*% c)
  se <- sqrt(V)
  
  # noncentrality parameter
  lambda <- true_theta / se
  
  # Satterthwaite df
  df_sat <- compute_satterthwaite_df(Sigma, c)
  
  # critical value from central t
  tcrit <- qt(1 - alpha/2, df = df_sat)
  
  # two-sided power using *noncentral* t
  power <- 
    pt(-tcrit, df = df_sat, ncp = lambda) + 
    (1 - pt(tcrit, df = df_sat, ncp = lambda))
  
  list(
    power  = power,
    df_sat = df_sat,
    se     = se,
    lambda = lambda
  )
}

estimate_power_replicated <- function(
    K = 50,
    B = 100,
    p_baseline = 0.40,
    betas,
    k_factor = 5,
    fit_interaction_order = 2, 
    icc_latent = 0.05,
    mean_m = 20,
    sd_m = 6.0553,
    alpha = 0.05,
    seed=123
) {
  set.seed(seed)
  # --------------------------------------------------------
  # STEP 1: Truth = FULL factorial model (all interactions)
  # --------------------------------------------------------
  test <- compute_factorial_probabilities(
    p_baseline       = p_baseline,
    betas            = betas,
    k_factor         = k_factor,
    interaction_order = k_factor
  )
  
  # random intercept variance
  logistic_var <- (pi^2)/3
  tau2 <- logistic_var * icc_latent / (1 - icc_latent)
  
  p_all_marginal <- compute_marginal_mean(
    eta  = test$eta_all,
    tau2 = tau2,
    method = "GH"
  )
  
  rho_marginal <- compute_marginal_icc(
    p          = p_all_marginal,
    icc_latent = icc_latent
  )
  
  W_res <- compute_W(
    p      = p_all_marginal,
    rho    = rho_marginal,
    mean_m = mean_m,
    sd_m   = sd_m
  )
  w_full <- W_res$w
  
  # --------------------------------------------------------
  # STEP 2: FULL X for variance model (matches manual block)
  # --------------------------------------------------------
  X_full <- cbind(1, test$design_full)
  colnames(X_full) <- c("(Intercept)", test$term_names)
  
  # Contrast lives in FULL X matrix
  c_vec <- rep(0, ncol(X_full))
  c_vec[2] <- 1     # component_1 always in column 2
  
  true_theta <- betas[1]
  
  # --------------------------------------------------------
  # STEP 3: Number of parameters in *fitted* model
  # --------------------------------------------------------
  p <- 1 + sum(choose(k_factor, 1:fit_interaction_order))
  
  # --------------------------------------------------------
  # STEP 4: Balanced expansion (matches manual block)
  # --------------------------------------------------------
  power_vec <- numeric(B)
  
  for (b in 1:B) {
    cl <- expand_design_for_clusters(X_full, w_full, K)
    X_cluster <- cl$X_cluster
    w_cluster <- cl$w_cluster
    
    Sigma <- compute_covariance_from_w(X_cluster, w_cluster)
    
    df <- max(K - p, 1)
    
    se_hat <- sqrt(t(c_vec) %*% Sigma %*% c_vec)
    lambda <- true_theta / se_hat
    tcrit  <- qt(1 - alpha/2, df)
    
    power_b <- pt(-tcrit, df = df, ncp = lambda) +
      (1 - pt(tcrit, df = df, ncp = lambda))
    
    power_vec[b] <- power_b
  }
  
  list(
    mean_power = mean(power_vec),
    sd_power   = sd(power_vec),
    all_power  = power_vec,
    df_used    = K - p,
    p          = p
  )
}

