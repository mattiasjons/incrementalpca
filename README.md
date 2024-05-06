# incrementalpca

R package to incrementally fit and update a PCA decomposition

# How to Install

Easiest way to install is:

```r
devtools::install_github("mattiasjons/incrementalpca", upgrade_dependencies = FALSE)
```

# How to use

A simple example for how to use the code follows below:

```r
library(incrementalpca)
d = 20
rho <- runif(d)
sigma <- d:1
df <- mvrnorm(100, rep(0, d), t(rho%*%diag(sigma))%*%t(rho))

inc.pca <- IncrementalDecomposition$new(df[1:20,], 10)
image(inc.pca$get_covariance(), ylim=c(1, 0))
image(inc.pca$get_precision(), ylim=c(1, 0))

inc.pca$partial_fit(df[21:50,], lambda = 0.98)
image(inc.pca$get_precision(), ylim=c(1, 0))
```
