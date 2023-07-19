# GVSSB
This is the repository for the paper A Generic Variational Spike-and-Slab Approach for Group Variable Selection.

To install the package, run the following code:

```
library(devtools)
install_github(repo="HowardGech/GVSSB")
```

If the package $\texttt{devtools}$ is not available in your R version,  download the zip file and unfold it into a user-specified folder (for example, GVSSB), and install it by:

```
install.packages("/path/to/folder/GVSSB", repos = NULL, type = "source")
```
The package includes three main functions:

GVSSB: This function computes the variational parameters for a given prior. It supports Gaussian prior, Laplace prior, and student T prior, and performs the variational Bayesian inference to estimate the model parameters.

predict.GVSSB: This function provides predictions of the response vector based on a fitted GVSSB model and a covariate matrix. It uses the estimated model parameters to generate predictions for new data points.

cv.GVSSB: This function performs cross-validation to select the prior with the smallest loss. It compares the performance of different priors by evaluating the loss function on validation data. Note that the function does not select the hyperparameters of the slab prior, as they are updated automatically within the GVSSB function.
