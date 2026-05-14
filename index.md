---
title: Power Estimation for Between-Cluster Randomized Factorial Trials with Binary Outcomes
---

# Overview

This page provides technical details for the method used in our poster:

> **Power Estimation for Between-Cluster Randomized Factorial Trials with Binary Outcomes**

We develop a fast approximation method for power calculation in
between-cluster randomized factorial trials (CRTs) with binary outcomes.

Key ideas:

- Distinguish **marginal** and **latent** ICCs and link them using an approximate relationship.
- Use a **cluster-level working linear model** for the logit of the cluster mean outcome.
- Approximate the variance of this logit either
  - analytically (delta method), or
  - via a fast Monte Carlo procedure.
- Plug these variances into a **GLS/Wald framework** to obtain power for main effects and interactions.

This page is intended as a technical appendix; the poster gives a short summary.

---

# 1. Model and ICC Definitions

## 1.1 Marginal ICC for Binary Outcomes

For two individuals \(i \ne j\) in the same cluster,

\[
\rho_{\text{marg}} = \mathrm{Corr}(Y_i, Y_j)
= \frac{\mathrm{Cov}(Y_i, Y_j)}
       {\sqrt{\mathrm{Var}(Y_i)\,\mathrm{Var}(Y_j)}}.
\]

For a binary outcome with marginal mean \(p\),

\[
\mathrm{Var}(Y_i) = p(1-p).
\]

\(\rho_{\text{marg}}\) is a probability-scale measure of within-cluster similarity and
automatically shrinks as \(p \to 0\) or \(p \to 1\).
This is the quantity that appears in most standard CRT power formulas.

## 1.2 Latent ICC in a Logistic Mixed Model

We work with a random-intercept logistic mixed model:

\[
  \text{logit}\{\Pr(Y_{ij}=1 \mid X_{ij}, u_j)\}
   = \beta_0 + \beta^\top X_{ij} + u_j,
  \qquad
  u_j \sim N(0,\tau^2).
\]

Equivalently, one can view a latent variable

\[
  Z_{ij} = X_{ij}^\top \beta + u_j + \varepsilon_{ij}, 
  \qquad
  \varepsilon_{ij} \sim \text{Logistic}(0,\pi^2/3),
\]

and define \(Y_{ij} = 1(Z_{ij} > 0)\).

The **latent ICC** is

\[
\rho_{\text{latent}}
  = \frac{\tau^2}{\tau^2 + \pi^2/3},
\]

which does not depend on prevalence and is natural for simulation studies that
specify \(\tau^2\) as a design parameter.

---

# 2. Working Model for Power Calculations

For power calculations, we approximate the GLMM by a **cluster-level
working linear model**:

\[
  Z_j = X_j \beta + \varepsilon_j, \qquad \varepsilon \sim (0, W),
\]

where \(Z_j\) is a transformed cluster-level summary (defined below),
\(X_j\) contains the factorial design indicators for cluster \(j\), and
\(W = \mathrm{diag}(w_1,\dots,w_J)\) is a diagonal working covariance (or weight)
matrix.

Under this GLS approximation, the covariance of \(\hat\beta\) is

\[
  \mathrm{Var}(\hat\beta)
  = (X^\top W^{-1} X)^{-1}.
\]

Power thus reduces to **choosing an appropriate summary \(Z_j\) and
approximating its variance** for each design cell / cluster.

---

# 3. The Key Idea: Logit of the Cluster Mean

We work at the cluster level and focus on the **logit of the cluster mean outcome**:

\[
  \bar Y_j = \frac{1}{m_j} \sum_{i=1}^{m_j} Y_{ij}, 
  \qquad
  Z_j = \text{logit}(\bar Y_j).
\]

We approximate the logistic mixed model by

\[
  \text{logit}(\bar Y_j)
  \;\approx\;
  X_j^\top \beta + \varepsilon_j,
\]

where \(\varepsilon_j\) captures the effects of clustering, the logistic
link, and unequal cluster sizes.

**Key reduction:** power estimation depends on approximating

\[
  \mathrm{Var}\{\text{logit}(\bar Y_j)\} = w_j
\]

for each cell \(j\). Once \(w_j\) is known, we assemble the working weight
matrix \(W\) and apply GLS/Wald theory.

We consider two complementary routes to approximate \(w_j\):

1. An **analytic (delta-method) approximation**, and  
2. A **Monte Carlo approximation** under the GLMM.

---

# 4. Analytic Approximation

## 4.1 Design and fixed effects

- Use an effect-coded \(2^k\) factorial design matrix \(X\) at the
  cluster level.
- Specify “true” fixed effects \(\beta\) (main effects and selected
  interactions) to reflect the underlying factorial structure.

For a given design cell with linear predictor \(\eta_j = \beta_0 + X_j^\top\beta\),
the conditional probability under the GLMM is

\[
  p_{\text{cond}, j}
  = \text{logit}^{-1}(\eta_j + u_j).
\]

## 4.2 Marginal mean and ICC

We marginalize over the random intercept:

\[
  p_j = \mathbb{E}_u \left[ \text{logit}^{-1}(\eta_j + u) \right],
  \qquad u \sim N(0,\tau^2),
\]

(evaluated numerically, if needed), and then convert the latent ICC to 
a **marginal ICC** using Goldstein’s approximation (see Section 6):

\[
  \rho_{\text{marg}, j}
  \approx 
  \frac{\tau^2\,p_j(1-p_j)}{1 + \tau^2\,p_j(1-p_j)}.
\]

## 4.3 Variance of the cluster mean

For a cluster of size \(m_j\) with marginal mean \(p_j\) and marginal ICC
\(\rho_{\text{marg}, j}\), the variance of the **cluster mean** is approximated by

\[
  \mathrm{Var}(\bar Y_j)
  \approx
  \frac{p_j(1-p_j)}{m_j}
  \left[1 + (m_j - 1)\rho_{\text{marg}, j}\right].
\]

If cluster sizes are unequal with mean \(\bar m\) and coefficient of
variation \(\text{CV}\), we may use an effective cluster size
\(m^\star = \bar m(1+\text{CV}^2)\) in place of \(m_j\).

## 4.4 Delta-method variance on the logit scale

Applying the delta method to the logit of the mean gives

\[
  \mathrm{Var}\{\text{logit}(\bar Y_j)\}
  \approx
  \frac{\mathrm{Var}(\bar Y_j)}{p_j^2(1-p_j)^2}
  =: w_j.
\]

We then form the **working weight matrix**

\[
  W = \mathrm{diag}(w_1,\dots,w_J),
\]

and obtain the working covariance of \(\hat\beta\) as

\[
  \Sigma = (X^\top W^{-1} X)^{-1}.
\]

---

# 5. Monte Carlo Approximation

Instead of relying on the delta method, we can approximate
\(\mathrm{Var}\{\text{logit}(\bar Y_j)\}\) **directly under the GLMM**.

For each factorial cell \(j\):

1. Draw a random intercept \(u_j \sim N(0,\tau^2)\).
2. Draw a cluster size \(m_j\) from a specified distribution.
3. Generate
   \[
     Y_{ij} \sim \text{Bernoulli}(\text{logit}^{-1}(\eta_j + u_j)),
     \quad i = 1,\dots,m_j.
   \]
4. Compute \(\bar Y_j\) and \(\text{logit}(\bar Y_j)\)
   (with a small continuity correction if \(\bar Y_j \in \{0,1\}\)).

Repeating this many times yields an empirical estimate

\[
  w_j = \mathrm{Var}\{\text{logit}(\bar Y_j)\}
\]

for each cell, which again defines the diagonal entries of \(W\).

Although this is more computationally intensive than the analytic route,
it is still very fast at the cluster level (under 1 second for typical
designs) and does not rely on normal approximations for \(\bar Y_j\).

---

# 6. Linking Latent and Marginal ICC

Power formulas for clustered binary data typically require the
**marginal ICC** \(\rho_{\text{marg}}\), while GLMM-based simulations
specify the **latent ICC** via \(\tau^2\).

Goldstein’s approximation provides a bridge between the two:

\[
  \rho_{\text{marg}}
  \approx
  \frac{\tau^2\,p(1-p)}{1 + \tau^2\,p(1-p)},
\]

where

\[
  \rho_{\text{latent}}
  = \frac{\tau^2}{\tau^2 + \pi^2/3}.
\]

This relationship shows how the marginal ICC shrinks toward zero as
\(p \to 0\) or \(p \to 1\) even when the latent ICC is fixed.

---

# 7. Power Calculation Given the Working Covariance

After expanding the \(2^k\) factorial design to \(K\) clusters using a
balanced allocation scheme, the GLS working covariance of \(\hat\beta\)
is

\[
  \Sigma = (X^\top W^{-1} X)^{-1}.
\]

For any contrast \(c^\top \beta\) (e.g., a main effect for one
component), we obtain the standard error

\[
  \mathrm{SE}(\hat\theta)
  = \sqrt{c^\top \Sigma\, c},
  \qquad \theta = c^\top \beta.
\]

We test

\[
  H_0 : \theta = 0 \quad \text{vs} \quad H_1 : \theta \neq 0
\]

using a Wald statistic

\[
  t = \frac{\hat\theta}{\mathrm{SE}(\hat\theta)},
\]

with a small-sample degrees-of-freedom correction

\[
  \text{df} = K - p,
\]

where \(p\) is the number of fixed-effect parameters in the model
(e.g., intercept, main effects, and 2-way interactions).

Under \(H_1\), power is computed from a noncentral t distribution with
noncentrality parameter

\[
  \lambda = \frac{\theta}{\mathrm{SE}(\hat\theta)}.
\]

---

# 8. Analytic vs. Monte Carlo Approximation

**Analytic approximation**

- Closed-form and very fast; well suited for design exploration.
- Relies on normal and delta-method approximations for \(\bar Y_j\) and
  \(\text{logit}(\bar Y_j)\).
- Accuracy improves with larger cluster sizes and moderate probabilities.

**Monte Carlo approximation**

- Directly reflects the GLMM data-generating process.
- Does not require normal approximations for cluster means.
- Slightly more computationally intensive, but still fast at the
  cluster-mean level.

**Practical recommendation**

- Use the analytic approximation as the default when average
  cluster sizes are moderate to large.
- Use the Monte Carlo approximation for small cluster sizes or as a
  sensitivity check.

---

# 9. Software and App

- **Interactive Shiny app (power calculator):**  
  \<link to your app\>

- **Code repository (this page):**  
  \<link to GitHub repo\>

---

# 10. Contact

For questions or comments, please contact:

- **Phuc Vu** – \<your email\>