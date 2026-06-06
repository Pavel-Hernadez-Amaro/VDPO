# VDPO 0.2.0

* Added support for partially observed functional data. This includes the new
  functions `ffpo()` and `ffpo_2d()` for defining partially observed functional
  terms in model formulae, `po_fit()` and `po_2d_fit()` for fitting the
  corresponding regression models, and the data generators
  `data_generator_po_1d()` and `data_generator_po_2d()`.
* `po_fit()` and `po_2d_fit()` now return the estimated functional coefficient
  together with its pointwise confidence intervals.
* Added `mfpca_vd()`, a multivariate functional principal component analysis for
  variable domain data, together with a `plot()` method for the estimated
  eigenfunctions and scores.
* Added vignettes for the partially observed models and for the variable domain
  multivariate functional principal component analysis.
* Completed the documentation and added runnable examples for the new functions.
* Fixed a bug in `plot_ci()` where the input object was referenced through the
  wrong variable, causing the function to error before producing a plot.

# VDPO 0.1.0

* Initial CRAN submission.
