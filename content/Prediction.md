---
title: Prediction and Inference
type: docs
prev: Model Fitting
math: true
---

## Sparse Group Linear Models
After model fitting, you can perform prediction using fitted `GVSSB` model and new data, simply by calling the `predict.GVSSB` function:

```r {filename="example - GVSSB prediction"}
n_test <- 50
X_test <- mvtnorm::rmvnorm(n_test, sigma=diag(p))
gvssb.predicted <- predict.GVSSB(gvssb.Laplace, X_test)
```
You may also get the volume of coverage region using the `coverage` function by specifying the fitted model, confidence level and inference type (`group` or `marginal`). If a coefficient vector is provided, it will also provide a vector of logical values indicating whether each coefficient (group) is selected as nonzero.

