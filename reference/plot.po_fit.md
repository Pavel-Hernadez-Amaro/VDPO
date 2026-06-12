# Plot method for partially observed functional regression fits

Displays an estimated functional coefficient of a `po_fit` object
together with its pointwise confidence band.

## Usage

``` r
# S3 method for class 'po_fit'
plot(x, beta_index = 1, ...)
```

## Arguments

- x:

  an object of class `po_fit`, as returned by
  [`po_fit`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_fit.md).

- beta_index:

  an integer selecting which functional coefficient to plot. Default is
  1.

- ...:

  currently ignored, included for compatibility with the generic.

## Value

A `ggplot2` object displaying the estimated coefficient and its
confidence band.

## See also

[`po_fit`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_fit.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(42)
  data <- data_generator_po_1d(N = 100)
  grid <- seq(0, 1, length.out = ncol(data$X_se))
  fit <- po_fit(y ~ ffpo(X_se, grid = grid), data = data)
  plot(fit)
}
#> Warning: unused arguments ignored for current version: N
#> Error in seq.default(0, 1, length.out = ncol(data$X_se)): 'length.out' must be of length 1
# }
```
