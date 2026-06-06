## Resubmission / update

This is an update of the VDPO package, which is already on CRAN (current
version 0.1.0). This submission updates it to version 0.2.0. It is not a new
package.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Summary of changes in this version

* Added support for partially observed functional data (`ffpo()`, `ffpo_2d()`,
  `po_fit()`, `po_2d_fit()` and the corresponding data generators).
* `po_fit()` and `po_2d_fit()` now return pointwise confidence intervals for
  the estimated functional coefficient.
* Added `mfpca_vd()`, a multivariate functional principal component analysis
  for variable domain data, with a `plot()` method.
* Added two new vignettes.

## Test environments

* Local: Windows 11, R 4.5.2.
* win-builder (devel).

## Downstream dependencies

There are currently no downstream dependencies for this package.
