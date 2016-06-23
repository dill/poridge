`poridge` - Principal co-Ordinate RIDGE regression
==================================================

An `R` package for principal coordinate ridge regression fitted by 'gam' (using `mgcv`). When the principal coordinates are defined by a relevant distance among functional predictors, this is a form of nonparametric scalar-on-function regression. Reiss et al. (2016) describe the approach and apply it to dynamic time warping distances among functional predictors.

An example (toy) analysis is shown [in the vignette](https://github.com/dill/poridge/blob/master/vignettes/PCoRR-vignette.pdf).

# Installation

The package can be installed via the `devtools` package directly from github using the following code:

```{r}
install.packages("devtools")
devtools::install_github("dill/poridge")
```

# Future development

The functions in `poridge` are available via the [`refund` package](https://github.com/refunders/refund)'s development branch. These methods will be available in the next [CRAN release](http://cran.r-project.org/package=refund) (c. July 2016).


# References

Reiss, P. T., Miller, D. L., Wu, P.-S., and Wen-Yu Hua, W.-Y. Penalized nonparametric scalar-on-function regression via principal coordinates. Under revision. [Preprint](https://works.bepress.com/phil_reiss/42/).

