# Generate high-signal 2D functional data with proper interface

Creates synthetic 2D functional data using the working high-signal
approach but with the same interface as data_generator_po_2d for
compatibility

## Usage

``` r
data_generator_high_signal(
  n = 100,
  grid_x = 20,
  grid_y = 20,
  intercept = 0.6,
  noise_sd = 0.25,
  response_type = "binomial",
  signal_strength = 2.5,
  n_missing = 0,
  min_distance_x = NULL,
  min_distance_y = NULL
)
```

## Arguments

- n:

  Number of samples to generate

- grid_x:

  Number of points in x-axis grid

- grid_y:

  Number of points in y-axis grid

- intercept:

  Target proportion for binomial or intercept for Gaussian

- noise_sd:

  Standard deviation of measurement noise

- response_type:

  Type of response ("binomial" or "gaussian")

- signal_strength:

  Strength of discriminative signal (default 2.5)

- n_missing:

  Not used (kept for compatibility)

- min_distance_x:

  Not used (kept for compatibility)

- min_distance_y:

  Not used (kept for compatibility)

## Value

A list containing simulated data with same structure as
data_generator_po_2d
