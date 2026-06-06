# Estimation of functional regression models for partially observed bidimensional functional data

The `po_2d_fit` function fits functional regression models for partially
observed bidimensional functional data, where each surface is only
observed over part of the common domain.

## Usage

``` r
po_2d_fit(formula, data, family = stats::gaussian(), offset = NULL)
```

## Arguments

- formula:

  a formula object with at least one `ffpo_2d` term.

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

An object of class `po_2d_fit`. It is a `list` containing the following
items:

- An item named `fit` of class `sop`. See
  [sop.fit](https://rdrr.io/pkg/SOP/man/sop.fit.html).

- An item named `Beta` which is a list with one entry per functional
  term, each containing the estimated coefficient surface (`beta`), its
  standard error (`se`), the lower and upper pointwise confidence limits
  (`lower`, `upper`) and the grids (`x`, `y`).

- An item named `intercept` which is the estimated intercept of the
  model.

- An item named `theta` which is the basis coefficient vector of the
  estimated bidimensional functional coefficient.

- An item named `covar_theta` which is the covariance matrix of the
  basis coefficients, used to build the pointwise confidence intervals.

- An item named `M` which holds the observed domain information for each
  functional term.

- An item named `ffpo_2d_evals` which is the result of the evaluations
  of the `ffpo_2d` terms in the formula.

## See also

[`ffpo_2d`](https://pavel-hernadez-amaro.github.io/VDPO/reference/ffpo_2d.md)

## Examples

``` r
# PARTIALLY OBSERVED BIDIMENSIONAL FUNCTIONAL DATA EXAMPLE
# \donttest{
# set seed for reproducibility
set.seed(123)

# generate example data with partially observed surfaces
sim  <- data_generator_po_2d(n = 30, grid_x = 8, grid_y = 8)
X    <- sim$noisy_surfaces_miss
y    <- sim$response
mp   <- sim$missing_points
mpts <- sim$miss_points

# Fit the model using an 'ffpo_2d' term for the partially observed surfaces.
# 'miss_points' and 'missing_points' describe the missing observations in the
# two complementary formats expected by 'ffpo_2d'.
fit  <- po_2d_fit(
  response ~ ffpo_2d(
    X = X, miss_points = mpts, missing_points = mp, nbasis = rep(5, 4)
  ),
  data = list(response = y, X = X)
)

# Inspect the structure of the returned object
str(fit, max.level = 1)
#> List of 7
#>  $ fit          :List of 15
#>   ..- attr(*, "class")= chr "sop"
#>  $ Beta         :List of 1
#>  $ intercept    : num 5.59
#>  $ theta        : num [1:25, 1] 2.78 -5.55 -13.88 -22.21 -30.54 ...
#>  $ covar_theta  : num [1:25, 1:25] 1191.1 622.7 54.2 -514.2 -1082.6 ...
#>  $ M            :List of 1
#>  $ ffpo_2d_evals:List of 1
#>  - attr(*, "class")= chr "po_2d_fit"
#>  - attr(*, "N")= int 30

# The basis coefficients of the functional coefficient can be accessed directly
fit$theta
#>              [,1]
#>  [1,]   2.7805011
#>  [2,]  -5.5501206
#>  [3,] -13.8806992
#>  [4,] -22.2111548
#>  [5,] -30.5415919
#>  [6,]   0.2185780
#>  [7,]  -4.0049825
#>  [8,]  -8.2285046
#>  [9,] -12.4519279
#> [10,] -16.6753384
#> [11,]  -2.3433450
#> [12,]  -2.4598444
#> [13,]  -2.5763101
#> [14,]  -2.6927011
#> [15,]  -2.8090850
#> [16,]  -4.9052681
#> [17,]  -0.9147063
#> [18,]   3.0758845
#> [19,]   7.0665258
#> [20,]  11.0571685
#> [21,]  -7.4671911
#> [22,]   0.6304318
#> [23,]   8.7280791
#> [24,]  16.8257526
#> [25,]  24.9234219

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
#> [1]   5.587516 -10.338622  -2.271412 -43.061656  38.527778
#> 
#> 
#> Estimated degrees of freedom:
#>       t_1       t_2 Total edf     Total 
#>         0         0         0         5 
#> 
#> R-sq.(adj) =  -0.284   Deviance explained = 11.4%  n = 30
#> 
#> Number of iterations: 1
# }
```
