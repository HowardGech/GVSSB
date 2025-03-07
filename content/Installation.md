---
title: Installation
type: docs
prev: /
next: Model Fitting
---

To install the GVSSB package, please ensure that the `devtools` package is installed in the R version you're using, then proceed with
```r {linenos=inline}
library(devtools)
install_github(repo="HowardGech/GVSSB")
```

If the package `devtools` is not available in your R version, you can also install it by source. Navigate to the [github](https://github.com/HowardGech/GVSSB/) page, download the zip file and unfold it into a user-specified folder (for example, GVSSB), and install it by:

```r
install.packages("/path/to/folder/GVSSB", repos = NULL, type = "source")
```
