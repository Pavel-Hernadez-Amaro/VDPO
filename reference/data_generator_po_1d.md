# Generate 1D functional data (current or legacy)

Provides the current 1D generator while keeping access to the previous
implementation via `version = "legacy"` for reproducibility.

## Usage

``` r
data_generator_po_1d(
  n = 100,
  grid_points = 100,
  noise_sd = 0.25,
  center = TRUE,
  rsq = 0.95,
  mu = 0.1,
  univariate = TRUE,
  response_type = c("gaussian", "binomial"),
  linear_predictor = c("rectangular", "trapezoidal", "linear"),
  n_missing = 1,
  min_distance = NULL,
  version = c("current", "legacy"),
  ...
)
```

## Arguments

- n:

  Number of samples to generate.

- grid_points:

  Number of points in the grid.

- noise_sd:

  Standard deviation of measurement noise.

- center:

  Whether to mean-center each curve.

- rsq:

  Desired R-squared value for the response.

- mu:

  Intercept term added to the linear predictor.

- univariate:

  If `TRUE`, generate a single functional predictor; otherwise generate
  two.

- response_type:

  Response distribution ("gaussian" or "binomial").

- linear_predictor:

  Integration approach for the linear predictor ("rectangular",
  "trapezoidal", or "linear").

- n_missing:

  Number of missing segments per curve.

- min_distance:

  Minimum length of missing segments (defaults to one fifth of the grid
  length).

- version:

  Choose `"current"` (default) or `"legacy"` implementation.

- ...:

  Additional arguments forwarded to the legacy implementation.

## Value

A list containing simulated curves (noisy and noiseless), missing point
indices, the coefficient functions, and the generated response.

## Examples

``` r
data <- data_generator_po_1d(n = 10)
data_legacy <- data_generator_po_1d(n = 10, version = "legacy", beta_type = "trig")
#> [1] 0.0008242013
```
