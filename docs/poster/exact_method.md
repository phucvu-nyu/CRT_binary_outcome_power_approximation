---
title: Method Overview
---

# Method overview (analytic and exact versions)

This page summarizes the main steps of our power approximation and explains how
the **analytic** and **exact** versions differ.

We assume the user specifies:

- A **baseline probability** \(p_{\text{baseline}}\),
- A set of **components/covariates** and their effects,
- A target **ICC** and cluster-size distribution,
- A total number of clusters \(K\).

---

## Step 1. Specify the factorial structure and effects

We describe the underlying effects in terms of a factorial structure:

- Choose the number of components \(k\) (e.g., 5 intervention components).
- Decide which main effects and interactions are nonzero.
- Encode these into a vector of coefficients (effect sizes).

Given:

- \(p_{\text{baseline}}\): the probability under a reference condition, and  
- the chosen effect sizes,

we can compute, for each combination of components (each “cell” of the design):

- A **linear predictor** \(\eta_j\),
- A corresponding **conditional probability**
  \(p_{\text{cond},j} = \text{logit}^{-1}(\eta_j)\),

under a logistic mixed model with random intercept.

At this stage, we have a set of probabilities for each cell, still conditional
on the random intercept being zero.

---

## Working data-generating model

Conceptually, we assume a **random-intercept logistic mixed model**:

\[
\text{logit}\{\Pr(Y_{ij} = 1)\}
= \beta_0 + \beta^\top X_{ij} + b_j,
\quad b_j \sim N(0,\tau^2),
\]

where

- \(Y_{ij}\) is the binary outcome for individual \(i\) in cluster \(j\),
- \(X_{ij}\) encodes the components/covariates for that individual/cluster,
- \(\beta_0, \beta\) are fixed effects chosen in Step 1,
- \(b_j\) is a cluster-specific random intercept with variance \(\tau^2\).

All subsequent steps (ICC mapping, marginal means, and variance approximations)
refer back to this working GLMM.

---

## Step 2. Connect ICC to the random-intercept variance

The user typically thinks in terms of an **ICC**, not a random-effect variance.
We work on the **latent** scale:

1. Start from a chosen **latent ICC** \(\rho_{\text{L}}\) (or convert a marginal
   ICC to \(\rho_{\text{L}}\) using Goldstein’s relationship).
2. Compute the corresponding random-intercept variance 
   \(
   \tau^2 = \frac{\pi^2}{3}\,\frac{\rho_{\text{L}}}{1 - \rho_{\text{L}}}.
   \)

This gives a fully specified random-intercept logistic model for each cell:
baseline, effects, and random-effect variance.

---

## Step 3. Obtain marginal means and marginal ICC (analytic setup)

To understand variability on the observed binary scale, we work with
**marginal means** and **marginal ICCs**:

1. Integrate out the random intercept to obtain **marginal means**
   \(
   p_{\text{marg},j}
   = \mathbb{E}\{\text{logit}^{-1}(\eta_j + b_j)\},
   \quad b_j \sim N(0,\tau^2),
   \)
   for each cell \(j\).
2. Use Goldstein’s approximation to convert the latent ICC \(\rho_{\text{L}}\)
   into a **marginal ICC** \(\rho_{\text{M},j}\) for each cell, given
   \(p_{\text{marg},j}\).

Now, for every design cell, we know:

- A marginal mean \(p_{\text{marg},j}\),
- A marginal ICC \(\rho_{\text{M},j}\).

---

## Step 4. Approximate the variance of the logit cluster mean

The central quantity for our power calculation is

\[
w_j \approx \mathrm{Var}\{\text{logit}(\bar Y_j)\}
\]

for each design cell \(j\), where \(\bar Y_j\) is the cluster mean outcome.

### Analytic version

In the **analytic** version we:

1. Summarize the cluster-size distribution by its mean \(\bar m\) and
   coefficient of variation, and convert this to an **effective cluster size**
   \(m^\star\).

2. For each cell \(j\), approximate the variance of the cluster mean
   \(\bar Y_j\) using a marginal-ICC CRT formula:
   \(\mathrm{Var}(\bar Y_j)
   \approx
   \frac{p_{\text{marg},j}(1 - p_{\text{marg},j})}{\bar m}
   \Big[ 1 + (m^\star - 1)\rho_{\text{M},j} \Big].
   \)
   

3. Apply the **delta method** to the logit function to get
   \(
   w_j
   \approx
   \frac{\mathrm{Var}(\bar Y_j)}
        {p_{\text{marg},j}^2 (1 - p_{\text{marg},j})^2}.
   \)

These \(w_j\) become the diagonal entries of a working weight matrix \(W\).

### Exact (GLMM-based) version

In the **exact** version, Step 4 is replaced by a dedicated Monte Carlo
approximation under the full random-intercept logistic model:

For each cell \(j\):

1. Draw a **cluster size** \(m_j\) from the chosen distribution
   (e.g. around the target mean).
2. Draw a **random intercept** \(b_j \sim N(0,\tau^2)\).
3. Use the GLMM to obtain a cluster-specific probability
   \(\Pr(Y_{ij}=1)\) for that cell.
4. Generate \(m_j\) binary outcomes in that cluster and compute the cluster
   mean \(\bar Y_j\).
5. Take \(\text{logit}(\bar Y_j)\) (with a small continuity correction).
6. Repeat many times and take the **empirical variance** of these logits; this
   is our estimate of \(w_j\).

Thus, both versions aim to approximate the same quantity
\(\mathrm{Var}\{\text{logit}(\bar Y_j)\}\), but they do so in different ways:
analytic formulas versus direct Monte Carlo under the GLMM.

---

## Step 5. From \(W\) to power (common to both)

Once we have \(w_j\) for each design cell, the remaining steps are identical for
the analytic and exact versions:

1. Form the diagonal matrix \(W = \mathrm{diag}(w_j)\).
2. Expand the set of cells to the desired **number of clusters** \(K\) using a
   balanced allocation of clusters to cells.
3. Use generalized least squares to obtain a **working covariance matrix**
   \(\Sigma\) for the estimated fixed effects.
4. For any contrast of interest (e.g. a main effect), compute
   \(\mathrm{SE}(\hat\theta)\) from \(\Sigma\) and obtain **power** from a
   noncentral \(t\) distribution with degrees of freedom \(K - p\), where
   \(p\) is the number of fixed-effect parameters in the fitted model.

The only difference between the analytic and exact methods is **how** \(w_j\) is
obtained in Step 4; the rest of the workflow is the same.