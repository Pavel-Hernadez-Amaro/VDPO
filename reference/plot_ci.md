# Plot Functional Curves with Confidence Intervals

Generates a plot of functional Beta estimates for specified curves,
along with their 95% confidence intervals. This function computes the
95% confidence intervals for each curve based on the covariance matrix
and the fitted values from the provided object. The resulting plot
includes estimated curves, confidence interval ribbons, and a legend
distinguishing the curves.

## Usage

``` r
plot_ci(object, beta_index = 1, curves)
```

## Arguments

- object:

  An object of class `'vd_fit'` or similar, containing the fitted model
  results, Beta estimates, and evaluation details.

- beta_index:

  An integer specifying which Beta coefficient matrix to use. Default is
  1.

- curves:

  A numeric vector specifying the indices of the curves (rows) to plot.

## Value

A `ggplot2` object displaying the Beta estimates and confidence
intervals for the specified curves.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  # Assuming `model_object` is an object of class 'vd_fit'
  plot_functional_curves_combined(model_object, beta = 1, curves = c(50, 70, 100))
}
} # }
```
