---
title: Model Fitting
type: docs
prev: Installation
next: Prediction
---

## Sparse Group Linear Models
The `GVSSB`{:.r} function aims at solving sparse group linear regression problems. For a given prior, including Gaussian prior, Laplace prior, and student T prior, it performs the variational Bayesian inference to estimate the model parameters. The main arguments include:
- `X`: A matrix of covariates in the linear regression problem, with rows being samples and columns being features.
- `Y`: An outcome vector or matrix with 1 column.
- `groups`: The group indicator vector with length the same length of features.
- `prior`: The slab part of the coefficient prior. To use the Cauchy distribution, set it to `"T"` and pass `nu=1'.
Other arguments can be found in its help function by calling `help(GVSSB)`{:.r}.