# Data generator function for the variable domain case

Generates a variable domain functional regression model

## Usage

``` r
data_generator_vd(
  N = 100,
  J = 100,
  nsims = 1,
  Rsq = 0.95,
  aligned = TRUE,
  multivariate = FALSE,
  beta_index = 1,
  use_x = FALSE,
  use_f = FALSE
)
```

## Arguments

- N:

  Number of subjects.

- J:

  Number of maximum observations per subject.

- nsims:

  Number of simulations per the simulation study.

- Rsq:

  Variance of the model.

- aligned:

  If the data that will be generated is aligned or not.

- multivariate:

  If TRUE, the data is generated with 2 functional variables.

- beta_index:

  Index for the beta.

- use_x:

  If the data is generated with x.

- use_f:

  If the data is generated with f.

## Value

A list containing the following components:

- y: `vector` of length N containing the response variable.

- X_s: `matrix` of non-noisy functional data for the first functional
  covariate.

- X_se: `matrix` of noisy functional data for the first functional
  covariate

- Y_s: `matrix` of non-noisy functional data for the second functional
  covariate (if multivariate).

- Y_se: `matrix` of noisy functional data for the second covariate (if
  multivariate).

- x1: `vector` of length N containing the non-functional covariate (if
  use_x is TRUE).

- x2: `vector` of length N containing the observed values of the smooth
  term (if use_f is TRUE).

- smooth_term: `vector` of length N containing a smooth term (if use_f
  is TRUE).

- Beta: `array` containing the true functional coefficients.

## Examples

``` r
# Basic usage with default parameters
sim_data <- data_generator_vd()

# Generate data with non-aligned domains
non_aligned_data <- data_generator_vd(N = 150, J = 120, aligned = FALSE)

# Generate multivariate functional data
multivariate_data <- data_generator_vd(N = 200, J = 100, multivariate = TRUE)

# Generate data with non-functional covariates and smooth term
complex_data <- data_generator_vd(
  N = 100,
  J = 150,
  use_x = TRUE,
  use_f = TRUE
)

# Generate data with a different beta function and R-squared value
custom_beta_data <- data_generator_vd(
  N = 80,
  J = 80,
  beta_index = 2,
  Rsq = 0.8
)

# Access components of the generated data
y <- sim_data$y # Response variable
X_s <- sim_data$X_s # Noise-free functional covariate
X_se <- sim_data$X_se # Noisy functional covariate
```
