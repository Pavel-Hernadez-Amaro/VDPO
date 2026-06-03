# Defining partially observed bidimensional functional data terms in VDPO formulae

Auxiliary function used to define `ffpo_2d` terms within `VDPO` model
formulae.

## Usage

``` r
ffpo_2d(X, miss_points, missing_points, nbasis = rep(15, 4), bdeg = rep(3, 4))
```

## Arguments

- X:

  partially observed bidimensional functional covariate `matrix`.

- miss_points, missing_points:

  `list` of missing observation points. See 'Details' for more
  information about the difference in structure between both.

- nbasis:

  number of basis to be used.

- bdeg:

  degree of the basis to be used.

## Value

The function is interpreted in the formula of a `VDPO` model. `list`
containing the following elements:

- `B_ffpo2d` design matrix.

- `Phi_ffpo2d` bidimensional B-spline basis used for the functional
  coefficient.

- `M_ffpo2d` the `missing_points` used as input in the function.

- `nbasis` number of the basis used.

## Details

The difference between miss_points and missing_points is the format in
which the data is presented.

`miss_points` is a `list` of `list`s where each inner list corresponds
to the observation points in the y-axis and contains the observation
points of the missing values for the x-axis. `miss_points` acts as a
guide for identifying and addressing missing observations in functional
data and is used for properly calculating the inner product matrix.

`missing_points` is a `list` where each element is a `matrix` containing
the missing observations points.
