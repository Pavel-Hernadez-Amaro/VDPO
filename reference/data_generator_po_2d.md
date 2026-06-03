# Generate 2D functional data for simulation studies

Creates synthetic 2D functional data with optional noise components and
different coefficient patterns. Uses Simpson's rule for accurate
integration.

## Usage

``` r
data_generator_po_2d(
  n = 20,
  grid_x = 20,
  grid_y = 20,
  noise_sd = 0.015,
  rsq = 0.95,
  intercept = 0.1,
  beta_type = c("saddle", "exp", "smooth", "sinusoidal", "peaks"),
  response_type = c("gaussian", "binomial"),
  linear_predictor = c("integral", "linear"),
  a1 = NULL,
  a2 = NULL,
  sub_response = 50,
  n_missing = 1,
  min_distance_x = NULL,
  min_distance_y = NULL
)
```

## Arguments

- n:

  Number of samples to generate.

- grid_x:

  Number of points in x-axis grid. Default is 20.

- grid_y:

  Number of points in y-axis grid. Default is 20.

- noise_sd:

  Standard deviation of measurement noise. Default is 0.015.

- rsq:

  Desired R-squared value for the response. Default is 0.95.

- intercept:

  Intercept added to the linear predictor (or target logit for binomial
  responses).

- beta_type:

  Type of coefficient surface ("saddle" or "exp"). Default is "saddle".

- response_type:

  Type of the response variable ("gaussian" or "binomial"). Default is
  "gaussian".

- linear_predictor:

  Integration approach for the linear predictor ("integral" or
  "linear").

- a1:

  Optional fixed value for first stochastic component. If provided, a2
  must also be provided.

- a2:

  Optional fixed value for second stochastic component. If provided, a1
  must also be provided.

- sub_response:

  Number of intervals for Simpson integration. Default is 50.

- n_missing:

  Number of holes in every curve.

- min_distance_x:

  Length of the holes in the x axis.

- min_distance_y:

  Length of the holes in the y axis.

## Value

A list containing:

- surfaces: List of n true (noiseless) surfaces

- noisy_surfaces: List of n observed (noisy) surfaces

- response: Vector of n response values

- grid_x: x-axis grid points

- grid_y: y-axis grid points

- beta: True coefficient surface

- stochastic_components: Matrix of a1 and a2 values used for each
  surface

## Examples

``` r
# Generate basic 2D functional data with default parameters
data <- data_generator_po_2d(n = 2)

# Generate data with custom grid size and Gaussian response
data <- data_generator_po_2d(n = 2, grid_x = 30, grid_y = 30, response_type = "gaussian")

# Generate data with binomial response and saddle-shaped coefficient surface
data <- data_generator_po_2d(n = 2, response_type = "binomial", beta_type = "saddle")
#> Warning: Maximum iterations reached without convergence. Final error: 0.15

# Generate data with fixed stochastic components
data <- data_generator_po_2d(n = 2, a1 = 0.1, a2 = -0.2)

# Introduce missing data with holes along curves
data <- data_generator_po_2d(n = 2, n_missing = 3, min_distance_x = 5, min_distance_y = 5)
#> Warning: numerical expression has 3 elements: only the first used
#> Warning: numerical expression has 3 elements: only the first used
#> Warning: numerical expression has 3 elements: only the first used
#> Warning: numerical expression has 3 elements: only the first used
```
