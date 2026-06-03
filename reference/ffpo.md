# Defining partially observed functional data terms in VDPO formulae

Auxiliary function used to define `ffpo` terms within `VDPO` model
formulae.

## Usage

``` r
ffpo(
  X,
  missing_points = NULL,
  grid,
  bidimensional_grid = FALSE,
  nbasis = c(30, 30),
  bdeg = c(3, 3),
  version = c("current", "legacy")
)
```

## Arguments

- X:

  partially observed functional covariate `matrix`.

- missing_points:

  observation points that were missing for each functional covariate
  `list`.

- grid:

  observation grid of the covariate.

- bidimensional_grid:

  boolean value that specifies if the grid should be treated as
  1-dimensional or 2-dimensional. The default value is `FALSE`
  (1-dimensional). See also 'Details'.

- nbasis:

  number of basis to be used.

- bdeg:

  degree of the basis to be used.

- version:

  Choose the `"current"` (default) or `"legacy"` implementation. The
  legacy version matches the previous `ffpo_old()` behavior.

## Value

the function is interpreted in the formula of a `VDPO` model. `list`
containing the following elements:

- `B_ffpo` design matrix.

- `Phi` B-spline basis used for the functional coefficient.

- `M` `vector` or `matrix` object indicating the observed domain of the
  data.

- `nbasis` number of the basis used.

## Details

When the same observation points are used for every functional
covariate, we end up with a vector observation grid. Imagine plotting
multiple curves, each representing a functional covariate, all measured
at the same time instances.

Conversely, if the observation points differ for each functional
covariate, we have a matrix observation grid. Picture a matrix where
each row represents a functional covariate, and the columns denote
distinct observation points. Varying observation points introduce
complexity, as each covariate might be sampled at different time
instances.

## See also

[`add_grid`](https://pavel-hernadez-amaro.github.io/VDPO/reference/add_grid.md)
