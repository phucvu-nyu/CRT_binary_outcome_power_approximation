---
title: Marginal vs Latent ICC
---

# Marginal vs latent ICC

## Marginal ICC

Defined on the observed binary scale:

\[
\rho_{\text{M}}
= \frac{\mathrm{Cov}(Y_{ij}, Y_{ik})}
       {\sqrt{\mathrm{Var}(Y_{ij}) \,\mathrm{Var}(Y_{ik})}}.
\]

- Depends on the outcome prevalence \(p\).
- Shrinks as \(p \to 0\) or \(p \to 1\).
- Common in design formulas for clustered binary data.

## Latent ICC

Random-intercept logistic mixed model:

\[
\text{logit}\{\Pr(Y_{ij}=1)\}
= \beta_0 + \beta^\top X_{ij} + b_j,
\quad b_j \sim N(0,\tau^2),
\]

with latent ICC

\[
\rho_{\text{L}}
= \frac{\tau^2}{\tau^2 + \pi^2/3}.
\]

- Does **not** depend on prevalence.
- Natural for specifying simulation models.

Our method links these two scales so that users can think in terms of either
marginal or latent ICC.