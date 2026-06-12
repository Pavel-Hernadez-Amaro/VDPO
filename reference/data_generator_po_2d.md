# Generate two-dimensional partially observed functional data

Simulates a scalar response together with partially observed functional
surfaces. The response is built from the integral of each surface
against a fixed coefficient surface, and a rectangular region of each
surface can be left unobserved.

## Usage

``` r
data_generator_po_2d(
  n = 100,
  grid_x = 20,
  grid_y = 20,
  intercept = 0.6,
  noise_sd = 0.25,
  response_type = c("binomial", "gaussian"),
  signal_strength = 2.5,
  beta_index = 1,
  n_missing = 0,
  min_distance_x = NULL,
  min_distance_y = NULL,
  verbose = FALSE
)
```

## Arguments

- n:

  Number of surfaces to generate.

- grid_x, grid_y:

  Number of grid points along each axis.

- intercept:

  Model intercept. For the binomial response it is used as the target
  proportion of successes.

- noise_sd:

  Standard deviation of the observation noise, relative to the standard
  deviation of each surface.

- response_type:

  Response distribution, either `"binomial"` (the default) or
  `"gaussian"`.

- signal_strength:

  Multiplier controlling the magnitude of the true coefficient surface.

- beta_index:

  Integer selecting the shape of the true coefficient surface, either
  `1` (the default, a wavy surface) or `2` (a simple inclined plane).

- n_missing:

  Number of unobserved rectangular regions per surface (default `0`,
  i.e. fully observed surfaces).

- min_distance_x, min_distance_y:

  Minimum size of the unobserved regions along each axis.

- verbose:

  If `TRUE`, print a short summary of the simulation. Defaults to
  `FALSE`.

## Value

A list with the true surfaces (`surfaces`), the noisy surfaces
(`noisy_surfaces`), the partially observed surfaces
(`noisy_surfaces_miss`) together with the missing point information
(`miss_points`, `missing_points`), the `response`, the true coefficient
surface (`beta`), the grids (`points_x`, `points_y`) and additional
simulation details.

## Examples

``` r
set.seed(123)
sim <- data_generator_po_2d(n = 20, grid_x = 10, grid_y = 10,
                            response_type = "gaussian")
str(sim, max.level = 1)
#> List of 16
#>  $ surfaces             :List of 20
#>  $ noisy_surfaces       :List of 20
#>  $ noisy_surfaces_miss  :List of 20
#>  $ miss_points          :List of 20
#>  $ missing_points       :List of 20
#>  $ response             : num [1:20] 0.666 0.467 0.5 0.784 0.723 ...
#>  $ intercept            : num 0.6
#>  $ points_x             : num [1:10] 0 0.111 0.222 0.333 0.444 ...
#>  $ points_y             : num [1:10] 0 0.111 0.222 0.333 0.444 ...
#>  $ beta                 : num [1:10, 1:10] 0.5 1.996 2.74 2.332 0.911 ...
#>  $ stochastic_components: num [1:20, 1:2] -0.441 -1.155 -1.152 1.268 1.111 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ functional_effects   : num [1:20] -0.0393 -0.1195 -0.1192 0.1531 0.1353 ...
#>  $ group_assignment     : num [1:20] 2 2 2 1 1 1 1 2 1 1 ...
#>  $ response_type        : chr "gaussian"
#>  $ signal_strength      : num 2.5
#>  $ diagnostics          :List of 2
```
