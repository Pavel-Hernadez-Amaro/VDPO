# Changelog

## VDPO 0.2.0

- Added support for partially observed functional data. This includes
  the new functions
  [`ffpo()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/ffpo.md)
  and
  [`ffpo_2d()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/ffpo_2d.md)
  for defining partially observed functional terms in model formulae,
  [`po_fit()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_fit.md)
  and
  [`po_2d_fit()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_2d_fit.md)
  for fitting the corresponding regression models, and the data
  generators
  [`data_generator_po_1d()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/data_generator_po_1d.md)
  and
  [`data_generator_po_2d()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/data_generator_po_2d.md).
- [`po_fit()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_fit.md)
  and
  [`po_2d_fit()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_2d_fit.md)
  now return the estimated functional coefficient together with its
  pointwise confidence intervals.
- Added
  [`mfpca_vd()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/mfpca_vd.md),
  a multivariate functional principal component analysis for variable
  domain data, together with a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
  the estimated eigenfunctions and scores.
- Added vignettes for the partially observed models and for the variable
  domain multivariate functional principal component analysis.
- Completed the documentation and added runnable examples for the new
  functions.
- Fixed a bug in
  [`plot_ci()`](https://pavel-hernadez-amaro.github.io/VDPO/reference/plot_ci.md)
  where the input object was referenced through the wrong variable,
  causing the function to error before producing a plot.

## VDPO 0.1.0

CRAN release: 2024-10-21

- Initial CRAN submission.
