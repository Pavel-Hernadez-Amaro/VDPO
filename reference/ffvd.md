# Defining variable domain functional data terms in vd_fit formulae

Auxiliary function used to define `ffvd` terms within `vd_fit` model
formulae. This term represents a functional predictor where each
function is observed over a domain of varying length. The formulation is
\\\frac{1}{T_i} \int \_1^{T_i} X_i(t)\beta(t,T_i)dt\\, where \\X_i(t)\\
is a functional covariate of length \\T_i\\, and \\\beta(t,T_i)\\ is an
unknown bivariate functional coefficient. The functional basis used to
model this term is the B-spline basis.

## Usage

``` r
ffvd(X, grid, nbasis = c(30, 50, 30), bdeg = c(3, 3, 3))
```

## Arguments

- X:

  variable domain functional covariate `matrix`.

- grid:

  observation points of the variable domain functional covariate. If not
  provided, it will be `1:ncol(X)`.

- nbasis:

  number of bspline basis to be used.

- bdeg:

  degree of the bspline basis used.

## Value

the function is interpreted in the formula of a `VDPO` model. `list`
containing the following elements:

- An item named `B` design matrix.

- An item named `X_hat` smoothed functional covariate.

- An item named `L_Phi` and `B_T` 1-dimensional marginal B-spline basis
  used for the functional coefficient.

- An item named `M` matrix object indicating the observed domain of the
  data.

- An item named `nbasis` number of basis used.

## Examples

``` r
# Generate sample data
set.seed(123)
data <- data_generator_vd(beta_index = 1, use_x = FALSE, use_f = FALSE)
X <- data$X_se

# Specifying a custom grid
custom_grid <- seq(0, 1, length.out = ncol(X))
ffvd_term_custom_grid <- ffvd(X, grid = custom_grid, nbasis = c(10, 10, 10))

# Customizing the number of basis functions
ffvd_term_custom_basis <- ffvd(X, nbasis = c(10, 10, 10))

# Customizing both basis functions and degrees
ffvd_term_custom <- ffvd(X, nbasis = c(10, 10, 10), bdeg = c(3, 3, 3))
```
