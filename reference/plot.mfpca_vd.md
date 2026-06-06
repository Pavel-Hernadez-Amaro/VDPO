# Plot method for variable domain multivariate FPCA

Displays either the estimated eigenfunctions of one variable at several
fixed domain lengths (superimposed lines) or the multivariate scores
colored by domain length.

## Usage

``` r
# S3 method for class 'mfpca_vd'
plot(
  x,
  type = c("eigenfunctions", "heatmap", "scores"),
  variable = 1,
  components = 1:2,
  domains = NULL,
  align_sign = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `mfpca_vd`.

- type:

  One of `"eigenfunctions"` (the default), `"heatmap"` or `"scores"`.
  `"eigenfunctions"` draws the chosen components at a few fixed domains
  as superimposed lines, `"heatmap"` shows one component across all
  domains, and `"scores"` plots the multivariate scores colored by
  domain length.

- variable:

  Index of the functional variable to display when
  `type = "eigenfunctions"` or `type = "heatmap"` (default `1`).

- components:

  Indices of the components to display (default `1:2`). For
  `type = "heatmap"` only the first one is used.

- domains:

  Domain lengths at which to draw the eigenfunctions. If `NULL`, a few
  values spanning the observed domains are used.

- align_sign:

  If `TRUE` (the default), the sign of each eigenfunction is aligned
  across domains so the curves overlay consistently.

- ...:

  Further graphical arguments.

## Value

Called for its side effect (a plot). Invisibly returns `NULL`.
