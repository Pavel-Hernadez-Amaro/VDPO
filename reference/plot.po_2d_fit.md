# Plot method for partially observed bidimensional functional regression fits

Displays an estimated bidimensional functional coefficient of a
`po_2d_fit` object as a heatmap over its two-dimensional domain.

## Usage

``` r
# S3 method for class 'po_2d_fit'
plot(x, beta_index = 1, ...)
```

## Arguments

- x:

  an object of class `po_2d_fit`, as returned by
  [`po_2d_fit`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_2d_fit.md).

- beta_index:

  an integer selecting which coefficient surface to plot. Default is 1.

- ...:

  currently ignored, included for compatibility with the generic.

## Value

A `ggplot2` object displaying the estimated coefficient surface.

## See also

[`po_2d_fit`](https://pavel-hernadez-amaro.github.io/VDPO/reference/po_2d_fit.md)
