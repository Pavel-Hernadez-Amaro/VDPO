# Multivariate functional principal component analysis for variable domain data

Performs a multivariate functional principal component analysis (MFPCA)
for functional data observed on subject-specific (variable) domains.
Each functional variable is first decomposed through a variable domain
FPCA, and the resulting univariate scores are combined into a
domain-varying multivariate decomposition.

## Usage

``` r
mfpca_vd(
  Data,
  Times = NULL,
  M_grid = NULL,
  Hz = 1,
  m_npcs = NULL,
  u_npcs = 5,
  k_m = 15,
  model_type = "gam"
)
```

## Arguments

- Data:

  A list (one entry per functional variable) of matrices with subjects
  in rows and observation points in columns, with `NA` on the unobserved
  part of each domain. All variables must have the same number of
  subjects. Currently two variables are supported.

- Times:

  A list of matrices, the same shape as `Data`, giving the actual
  observation time of each measurement. Mutually exclusive with
  `M_grid`.

- M_grid:

  A numeric vector of length \\N\\ with the domain length of each
  subject, used when the observations are equidistant. Mutually
  exclusive with `Times`.

- Hz:

  Sampling rate used to build the equidistant grid when `M_grid` is
  supplied (default `1`).

- m_npcs:

  Number of multivariate principal components to retain. If `NULL`, the
  number explaining more than 90 percent of the variance is used.

- u_npcs:

  Number of univariate principal components retained per variable
  (default `5`).

- k_m:

  Dimension of the basis used to model the score covariances along the
  domain (default `15`).

- model_type:

  Either `"gam"` (the default) or `"sop"`, the method used to smooth the
  score covariances along the domain.

## Value

A list with the following items:

- `scores_m`:

  Matrix of multivariate scores.

- `efunctions_m`:

  Multivariate eigenfunctions, as a list over domains, each a list over
  variables.

- `efunctions_u`:

  Univariate eigenfunctions, as a list over variables.

- `scores_u`:

  Matrix of univariate scores.

- `evalues_u`:

  Univariate eigenvalues, as a list over variables.

- `evalues_m`:

  Multivariate eigenvalues, as a list over domains.

- `var_u`:

  Cumulative variance explained by the univariate components.

- `mean_model`:

  The fitted univariate mean models (one per variable).

- `M_grid`:

  The domain grid used.

- `argvals_u`:

  The observation times of each subject, as a list over variables, each
  a list over subjects.

## Details

The observation times can be supplied in two mutually exclusive ways.
With `Times`, each variable carries a matrix (subjects in rows,
observation points in columns) holding the actual observation time of
every measurement, and the domain of subject \\i\\ is taken as the
maximum observed time. With `M_grid`, the observations are assumed
equidistant and `M_grid` is a vector of length \\N\\ giving the domain
length of each subject, so subject \\i\\ is evaluated on
`seq(0, M_grid[i], by = 1 / Hz)`. Exactly one of `Times` or `M_grid`
must be provided.
