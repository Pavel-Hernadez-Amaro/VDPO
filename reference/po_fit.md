# Estimation of functional regression models for partially observed functional data

The `po_fit` function fits functional regression models for partially
observed functional data, where each curve is only observed over part of
the common domain.

## Usage

``` r
po_fit(formula, data, family = stats::gaussian(), offset = NULL)
```

## Arguments

- formula:

  a formula object with at least one `ffpo` term.

- data:

  a `list` object containing the response variable and the covariates as
  the components of the list.

- family:

  a `family` object specifying the distribution from which the data
  originates. The default distribution is
  [`gaussian`](https://rdrr.io/r/stats/family.html).

- offset:

  an offset vector. The default value is `NULL`.

## Value

An object of class `po_fit`. It is a `list` containing the following
items:

- An item named `fit` of class `sop`. See
  [sop.fit](https://rdrr.io/pkg/SOP/man/sop.fit.html).

- An item named `Beta` which is a list with one `data.frame` per
  functional term, each containing the grid (`t`), the estimated
  functional coefficient (`beta`), its standard error (`se`) and the
  lower and upper limits of the pointwise confidence interval (`lower`,
  `upper`).

- An item named `intercept` which is the estimated intercept of the
  model.

- An item named `theta` which is the basis coefficient vector of the
  estimated functional coefficient.

- An item named `covar_theta` which is the covariance matrix of the
  basis coefficients, used to build the pointwise confidence intervals.

- An item named `M` which holds the observed domain information for each
  functional term.

- An item named `ffpo_evals` which is the result of the evaluations of
  the `ffpo` terms in the formula.

## See also

[`ffpo`](https://pavel-hernadez-amaro.github.io/VDPO/reference/ffpo.md)

## Examples

``` r
# PARTIALLY OBSERVED FUNCTIONAL DATA EXAMPLE

# set seed for reproducibility
set.seed(123)

# generate example data with partially observed curves
sim  <- data_generator_po_1d(n = 50, grid_points = 60)
X    <- sim$noisy_curves_miss
y    <- sim$response
grid <- sim$grid
mp   <- sim$missing_points

# Fit the model using an 'ffpo' term for the partially observed covariate.
# 'nbasis' sets the number of basis functions for the data reconstruction
# and for the functional coefficient, respectively.
fit  <- po_fit(
  response ~ ffpo(X = X, missing_points = mp, grid = grid, nbasis = c(20, 20)),
  data = list(response = y, X = X, grid = grid, missing_points = mp)
)

# Inspect the structure of the returned object
str(fit, max.level = 1)
#> List of 7
#>  $ fit        :List of 15
#>   ..- attr(*, "class")= chr "sop"
#>  $ Beta       :List of 1
#>  $ intercept  : num 0.123
#>  $ theta      : num [1:20, 1] -0.03136 -0.01256 0.00602 0.02258 0.03373 ...
#>  $ covar_theta: num [1:20, 1:20] 8.72e-04 4.35e-04 1.23e-04 -1.74e-05 -3.27e-05 ...
#>  $ M          :List of 1
#>  $ ffpo_evals :List of 1
#>  - attr(*, "class")= chr "po_fit"
#>  - attr(*, "N")= int 50

# The estimated intercept and the basis coefficients can be accessed directly
fit$intercept
#> [1] 0.1229534
fit$theta
#>               [,1]
#>  [1,] -0.031363841
#>  [2,] -0.012559819
#>  [3,]  0.006023635
#>  [4,]  0.022584795
#>  [5,]  0.033734844
#>  [6,]  0.036421863
#>  [7,]  0.033038168
#>  [8,]  0.026912137
#>  [9,]  0.015454328
#> [10,] -0.002773189
#> [11,] -0.023712627
#> [12,] -0.042205336
#> [13,] -0.057212897
#> [14,] -0.068095604
#> [15,] -0.073799895
#> [16,] -0.072385012
#> [17,] -0.062012250
#> [18,] -0.044650855
#> [19,] -0.025073615
#> [20,] -0.005355858

# A summary of the underlying fit can be obtained using the summary function
summary(fit)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> 
#> Formula:
#> NULL
#> 
#> 
#> Fixed terms: 
#> [1]  0.1229534 -0.0116623  0.1245468
#> 
#> 
#> Estimated degrees of freedom:
#> Total edf     Total      <NA> 
#>    3.6034    3.6034    6.6034 
#> 
#> R-sq.(adj) =  0.732   Deviance explained = 81.9%  n = 50
#> 
#> Number of iterations: 1
```
