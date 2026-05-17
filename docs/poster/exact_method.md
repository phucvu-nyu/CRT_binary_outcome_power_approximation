---
title: Exact Method Implementation
---

# Exact GLMM-based method

We implement a **Monte Carlo approximation** to
\(\mathrm{Var}\{\text{logit}(\bar Y_j)\}\) under the full logistic mixed model.

Key steps (see `new_method_exact.R`):

1. **Factorial truth**  
   - Build a \(2^k\) factorial design (effect coded \(\pm 1\)).
   - Specify main effects and interactions via `betas`.
   - Calibrate the intercept so that the all-\(-1\) cell has baseline
     probability \(p_{\text{baseline}}\).

2. **Random intercept variance**  
   - Specify a latent ICC \(\rho_{\text{L}}\).  
   - Convert to variance
     \(\tau^2 = \frac{\pi^2}{3}\,\rho_{\text{L}}/(1-\rho_{\text{L}})\).

3. **Variance of logit cluster means**  
   - For each factorial cell \(j\), repeatedly:
     - Draw a random intercept \(b_j \sim N(0,\tau^2)\).
     - Draw a cluster size \(m_j\) from a chosen distribution.
     - Generate \(Y_{ij} \sim \text{Bernoulli}(\text{logit}^{-1}(\eta_j + b_j))\).
     - Compute \(\bar Y_j\) and \(\text{logit}(\bar Y_j)\).
   - Use the empirical variance of these logits as \(w_j\).

4. **Power**  
   - Build \(W = \text{diag}(w_1,\dots,w_J)\).  
   - Expand the factorial design to \(K\) clusters.  
   - Form \(\Sigma = (X^\top W^{-1} X)^{-1}\).  
   - Compute power for a chosen contrast (e.g., main effect for component 1)
     using a noncentral \(t\) with df \(K - p\).

All of this logic is wrapped in `estimate_power_replicated_glmm_mc()` in
`new_method_exact.R`.