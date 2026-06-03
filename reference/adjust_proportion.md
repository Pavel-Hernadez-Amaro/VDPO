# Iteratively adjust intercept to achieve target proportion in binomial simulation

This function uses an iterative approach to find the appropriate
intercept value that produces a desired proportion of 1s in binomial
response simulation. It works by adjusting the intercept on the log-odds
scale using adaptive damping to prevent overshooting due to the
nonlinear logistic transformation.

## Usage

``` r
adjust_proportion(
  target_prop,
  functional_effects,
  max_iter = 15,
  tolerance = 0.03,
  verbose = FALSE
)
```

## Arguments

- target_prop:

  Desired proportion of 1s in the response. Must be between 0 and 1
  (exclusive).

- functional_effects:

  Numeric vector of functional effects (e.g., from 2D integration of
  surfaces). These represent the variability around the baseline
  intercept.

- max_iter:

  Maximum number of iterations for adjustment. Default is 15.

- tolerance:

  Convergence tolerance for the difference between target and achieved
  proportion. Default is 0.03.

- verbose:

  Logical indicating whether to print iteration progress. Default is
  FALSE.

## Value

A list containing:

- y: Binary response vector of length equal to `functional_effects`

- intercept: Final adjusted intercept value

- final_prop: Achieved proportion of 1s in the response

- iterations: Number of iterations used

- converged: Logical indicating whether convergence was achieved

- final_error: Final absolute error between target and achieved
  proportion

## Details

The function works by:

1.  Starting with an initial intercept based on `qlogis(target_prop)`

2.  Computing probabilities using
    `plogis(intercept + functional_effects)`

3.  Generating binary outcomes using
    [`rbinom()`](https://rdrr.io/r/stats/Binomial.html)

4.  Adjusting the intercept based on the error between target and
    achieved proportions

5.  Using adaptive damping (0.3 for large errors, 0.5 for medium, 0.7
    for small)

The adjustment formula is:
`adjustment = (qlogis(target_prop) - qlogis(achieved_prop)) * damping_factor`

This approach works in log-odds space to prevent probabilities from
exceeding the unit interval and provides robust control over response
proportions in functional regression simulation studies.

## See also

[`data_generator_po_2d`](https://pavel-hernadez-amaro.github.io/VDPO/reference/data_generator_po_2d.md)
for using this function in 2D functional data simulation.

## Examples

``` r
# Basic usage with simulated functional effects
set.seed(123)
effects <- rnorm(100, mean = 0, sd = 0.5)
result <- adjust_proportion(target_prop = 0.3, functional_effects = effects)
cat("Achieved proportion:", result$final_prop, "\n")
#> Achieved proportion: 0.27 
cat("Converged in", result$iterations, "iterations\n")
#> Converged in 1 iterations

# Usage with 2D functional regression effects
# Assuming you have computed 2D integral effects from surfaces
# integral_effects <- compute_2d_integrals(surfaces, beta_surface)
# result <- adjust_proportion(0.4, integral_effects, tolerance = 0.01)

# Check for convergence issues
if (!result$converged) {
  warning("Adjustment did not converge. Final error: ", result$final_error)
}
```
