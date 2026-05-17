---
title: General Idea
---

# General idea of the method

Our power calculation is built on a **working linear model** at the **cluster
level**, not on the individual binary outcomes directly.

---

## 1. Work with cluster means, not individuals

For each factorial cell \(j\) (a particular combination of components), we look
at the **cluster mean outcome**

\[
\bar Y_j = \frac{1}{m_j} \sum_{i=1}^{m_j} Y_{ij},
\]

and then take its logit:

\[
Z_j = \text{logit}(\bar Y_j).
\]

The key approximation is that these logits follow a linear model:

\[
Z_j \approx X_j \beta + \varepsilon_j,
\]

where

- \(X_j\) encodes the factorial components (main effects + interactions),
- \(\beta\) are the fixed effects of interest,
- \(\varepsilon_j\) collects **all** of the complications:
  clustering, nonlinearity of the logit link, unequal cluster sizes, etc.

We do **not** require \(Z_j\) to be exactly normal; we only need a good handle
on its variance.

---

## 2. Reduce the problem to a variance

Once we accept the working model

\[
Z_j \approx X_j \beta + \varepsilon_j,
\]

standard GLS/Wald theory tells us that

\[
\mathrm{Var}(\hat\beta) \approx (X^\top W^{-1} X)^{-1},
\]

where \(W\) is a diagonal matrix with entries

\[
W = \mathrm{diag}(w_1,\dots,w_J),
\quad
w_j \approx \mathrm{Var}\{Z_j\}
      = \mathrm{Var}\{\text{logit}(\bar Y_j)\}.
\]

So **all of the power calculation** boils down to:

> How do we approximate  
> \(\mathrm{Var}\{\text{logit}(\bar Y_j)\}\) for each factorial cell \(j\)?

If we can approximate these variances well, we get:

1. \(W = \mathrm{diag}(w_j)\),
2. \(\Sigma = (X^\top W^{-1} X)^{-1}\),
3.  For any contrast \(c^\top \beta\) (e.g., a main effect),
   \(\mathrm{SE}(\hat\theta) \approx \sqrt{c^\top \Sigma c}\),
4. And then **power** from a noncentral \(t\) test with df \(K - p\).

---

## 3. Two routes for the variance

In the rest of the site we describe two ways to approximate
\(\mathrm{Var}\{\text{logit}(\bar Y_j)\}\):

1. An **analytic** approximation using marginal means and ICCs.
2. An **exact GLMM-based Monte Carlo** approximation that works directly under
   the random-intercept logistic model.

The poster and this website focus on the second route: the exact Monte Carlo
approximation implemented in `new_method_exact.R`.