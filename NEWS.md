# VDPO 0.2.0

* Added support for partially observed functional data. This includes the new
  functions `ffpo()` and `ffpo_2d()` for defining partially observed functional
  terms in model formulae, `po_fit()` and `po_2d_fit()` for fitting the
  corresponding regression models, and the data generators
  `data_generator_po_1d()` and `data_generator_po_2d()`.
* Completed the documentation and added runnable examples for the partially
  observed functions.
* Fixed a bug in `plot_ci()` where the input object was referenced through the
  wrong variable, causing the function to error before producing a plot.

# VDPO 0.1.0

* Initial CRAN submission.
