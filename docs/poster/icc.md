---
title: Marginal vs Latent ICC
---

# Marginal vs latent ICC

Cluster trials with binary outcomes often talk about “the ICC,” but there are
**two different ICCs** in play. They answer slightly different questions. This page is dedicated to clarifying the two ICC definitions used in clustered
binary trials—**marginal** and **latent**—and to explaining how we link them
using **Goldstein’s approximation**. This link is what lets us move cleanly
between design-scale quantities and the GLMM used for simulation and power.

---

## 1. Marginal ICC: correlation on the observed scale

Imagine two individuals \(i \neq k\) in the same cluster \(j\).  
The **marginal ICC** is the correlation between their binary outcomes:

\[
\rho_{\text{M}}
= \frac{\mathrm{Cov}(Y_{ij}, Y_{kj})}
       {\sqrt{\mathrm{Var}(Y_{ij}) \,\mathrm{Var}(Y_{kj})}} .
\]

For binary outcomes, \(\mathrm{Var}(Y_{ij}) = p_j(1-p_j)\), so:

- \(\rho_{\text{M}}\) **depends on the prevalence** \(p_j\).
- As \(p_j \to 0\) or \(p_j \to 1\), the variance \(p_j(1-p_j)\) shrinks and so
  does \(\rho_{\text{M}}\).

**Interpretation**

- “If one person in the cluster has \(Y=1\), how much more likely is another
  person in the same cluster to have \(Y=1\), on the *observed* 0/1 scale?”
- This is the ICC that appears in many textbook sample size formulas for
  clustered binary data.

---

## 2. Latent ICC: correlation on an underlying continuous scale

For modeling and simulation, we often use a random–intercept logistic mixed
model:

\[
\text{logit}\{\Pr(Y_{ij}=1)\}
= \beta_0 + \beta^\top X_{ij} + b_j,
\quad b_j \sim N(0,\tau^2).
\]

Equivalently, we can think of a continuous **latent variable**

\[
Z_{ij} = X_{ij}\beta + b_j + \varepsilon_{ij},
\quad \varepsilon_{ij} \sim \text{Logistic}(0,\pi^2/3),
\]

with \(Y_{ij} = 1\) if \(Z_{ij} > 0\).

On this latent scale, the **latent ICC** is

\[
\rho_{\text{L}}
= \frac{\tau^2}{\tau^2 + \pi^2/3}.
\]

Key features:

- \(\rho_{\text{L}}\) is **constant across prevalences**.
- It is determined entirely by the random–effect variance \(\tau^2\).
- It is very convenient for simulations: you can choose \(\rho_{\text{L}}\),
  compute \(\tau^2\), and plug it into the GLMM.

---

## 3. Why we need both

- Trial **design formulas** are written in terms of **marginal ICC**
  \(\rho_{\text{M}}\) on the observed scale.
- **Simulation models** and mixed‑model software work naturally with the
  **latent ICC** \(\rho_{\text{L}}\).

If we are not careful, we can:

- Specify a trial assuming a certain marginal ICC,
- But simulate or analyze it using a mismatched latent ICC.

---

## 4. Linking the two: Goldstein’s approximation

Let \(p\) be the marginal mean outcome in a given cell, and let \(\tau^2\) be the
random–intercept variance from the logistic mixed model.

Goldstein’s approximation gives a direct relationship between the **latent** and
**marginal** ICCs:

\[
\rho_{\text{M}}
\approx
\frac{\tau^2\,p(1-p)}{1 + \tau^2\,p(1-p)},
\]

where the latent ICC is

\[
\rho_{\text{L}}
= \frac{\tau^2}{\tau^2 + \pi^2/3}.
\]

This approximation has two key implications:

- For fixed \(\tau^2\) (fixed \(\rho_{\text{L}}\)), the marginal ICC
  \(\rho_{\text{M}}\) decreases as \(p \to 0\) or \(p \to 1\).
- Given either \(\rho_{\text{L}}\) or \(\rho_{\text{M}}\) (plus \(p\)), we can
  **translate between the two scales**.

