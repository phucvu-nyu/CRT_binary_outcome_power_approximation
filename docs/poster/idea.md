---
title: General Idea
---

# General idea of the method

We work at the **cluster level** and focus on the **logit of the cluster mean**.

- Start from a logistic mixed model (with random intercept).
- For each factorial cell, define the **cluster mean** \(\bar Y_j\).
- Approximate
  \[
    \text{logit}(\bar Y_j) \approx X_j \beta + \varepsilon_j,
  \]
  where \(\varepsilon_j\) captures clustering, nonlinearity, and unequal sizes.

Given an approximation to
\(\mathrm{Var}\{\text{logit}(\bar Y_j)\}\) for each cell, we can:

1. Form a diagonal weight matrix \(W\).
2. Compute \(\Sigma = (X^\top W^{-1} X)^{-1}\).
3. Use standard **Wald / noncentral \(t\)** theory for a contrast \(c^\top \beta\)
   to get power.

The only hard part is approximating
\(\mathrm{Var}\{\text{logit}(\bar Y_j)\}\) well for the binary clustered case.