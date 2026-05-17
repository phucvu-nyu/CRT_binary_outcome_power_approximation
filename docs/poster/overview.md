---
title: Overview and Challenges
---

# Power estimation for factorial cluster trials

We consider power estimation for **between-cluster randomized factorial trials
with binary outcomes**.

Existing tools for factorial trials mostly target **continuous** outcomes and
often require programming. For binary outcomes in clustered designs, power is
typically obtained by **full simulation**, which can be:

- Computationally intensive
- Sensitive to poorly chosen ICCs
- Prone to inflated type I error without small-sample corrections

Our goal is to provide a **fast, simulation-accurate** power method that:

- Works for multi-component factorial designs
- Handles binary clustered outcomes
- Is simple enough to embed in a user-facing app