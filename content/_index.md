---
title: ""
# toc: false
type: docs
math: true
next: Installation
---
![](GVSSB_icon.svg)

**GVSSB** a parameter-expanded Coordinate-Ascent Variational Inference algorithm for high-dimensional linear regression models with grouped variables under spike-and-slab priors in R language. The sparse linear models can be formulated as
$$\bm{Y} = \sum_{j = 1}^G\bm{X}_j\bm{\beta}+\bm{\epsilon},\quad \bm{\epsilon}\sim N(\bm{0}, \sigma \bm{I}_n)$$
where the coefficient $\bm{\beta}$ is sparse with most of its elements being zero. It is also capable of solving the sparse additive model
$$Y_i=b+\sum_{j=1}^pf_j(X_{ij})+\epsilon_i,\quad \epsilon_i\sim N(0,\sigma_\star^2)$$
with most of $f_j$'s being zero functions by spline fitting. The smothness of splines are regularized to prevent overfitting. For more information please refer [here](https://arxiv.org/abs/2309.16855).

<!-- {{< cards >}}
  {{< card link="Documentation" title="Docs" icon="book-open" >}}
{{< /cards >}}
{{< cards >}}
  {{< card link="s" title="Docs" icon="book-open" >}}
{{< /cards >}} -->
