#' Add missing values to curves
#'
#' @param X List of curves
#' @param n_missing Number of holes in every curve
#' @param min_distance Length of the holes
#'
#' @return List containing curves with missing values and missing points information
#'
#' @noRd
add_miss1d <- function(X, n_missing = 1, min_distance = 5) {
  N <- length(X)
  missing_points <- miss_points <- vector(mode = "list", length = N)

  n_points <- length(X[[1]])

  for (i in 2:N) {
    if (length(X[[i]]) != n_points) {
      stop("The length of all curves must be the same.", call. = FALSE)
    }
  }

  for (j in seq_along(missing_points)) {
    if (n_missing >= 1) {
      x_missing <- sort(sample(1:(n_points - min_distance), n_missing))

      missing_ranges <- lapply(x_missing, function(x) x:(x + min_distance - 1))
      missing_points[[j]] <- unlist(missing_ranges)

      X[[j]][missing_points[[j]]] <- NA
    }

    miss_points[[j]] <- which(is.na(X[[j]]))
  }

  list(
    X_miss = X,
    miss_points = miss_points,
    missing_points = missing_points
  )
}

#' Generate 1D functional data for simulation studies
#'
#' Creates synthetic 1D functional data with optional noise components and different
#' coefficient patterns. Uses trapezoidal rule for integration.
#'
#' @param n Number of samples to generate
#' @param grid_points Number of points in the grid. Default is 100
#' @param noise_sd Standard deviation of measurement noise. Default is 0.015
#' @param rsq Desired R-squared value for the response. Default is 0.95
#' @param beta_type Type of coefficient function ("sin" or "gaussian"). Default is "sin"
#' @param n_missing Number of missing segments per curve. Default is 1
#' @param min_distance Minimum length of missing segments. Default is NULL (auto-calculated)
#'
#' @return A list containing:
#' \itemize{
#'   \item curves: List of n true (noiseless) curves
#'   \item noisy_curves: List of n observed (noisy) curves
#'   \item noisy_curves_miss: List containing curves with missing values
#'   \item response: Vector of n response values
#'   \item grid: Grid points
#'   \item beta: True coefficient function
#'   \item stochastic_components: Vector of a values used for each curve
#' }
#'
generate_1d_po_functional_data <- function(
    n = 100,
    grid_points = 100,
    noise_sd = 0.015,
    rsq = 0.95,
    beta_type = c("sin", "gaussian"),
    n_missing = 1,
    min_distance = NULL
) {

  beta_type <- match.arg(beta_type)

  if (is.null(min_distance)) {
    min_distance <- round(1/5 * grid_points)
  }

  # Create grid points
  t <- seq(0, 1, length.out = grid_points)

  # Initialize storage
  curves <- vector("list", n)
  noisy_curves <- vector("list", n)
  stochastic_components <- vector("list", n)
  nu <- numeric(n)

  # Helper function to generate curve with given parameters
  generate_curve <- function(t, a1, a2, b1, noise_sd) {
    true_curve <- a1 * sin(2 * pi * t) +
      a2 * cos(4 * pi * t) +
      b1 * t +
      exp(-t) +
      1

    list(
      curve_true = true_curve,
      curve_noisy = true_curve + rnorm(length(t), 0, noise_sd)
    )
  }


  # Generate coefficient function
  beta <- if (beta_type == "sin") {
    sin(2 * pi * t) * cos(pi * t)
  } else {
    dnorm(t, mean = 0.5, sd = 0.15) / dnorm(0.5, mean = 0.5, sd = 0.15)
  }

  # Generate data and compute response
  for (i in 1:n) {
    # Use fixed value if provided, otherwise generate random component
    stochastic_components[[i]] <- rnorm(3, 0, 0.2)

    # Generate curve
    curve_data <- generate_curve(
      t,
      stochastic_components[[i]][1],
      stochastic_components[[i]][2],
      stochastic_components[[i]][3],
      noise_sd
    )

    curves[[i]] <- curve_data$curve_true
    noisy_curves[[i]] <- curve_data$curve_noisy

    # Compute integral using trapezoidal rule
    integrand <- curves[[i]] * beta
    nu[i] <- 0.5 * sum(diff(t) * (integrand[-1] + integrand[-length(integrand)]))
  }

  # Generate response with desired R-squared
  var_e <- (1/rsq - 1) * var(nu)
  response <- nu + rnorm(n, 0, sqrt(var_e))

  # Add missing values
  noisy_curves_miss <- add_miss1d(
    noisy_curves,
    n_missing = n_missing,
    min_distance = min_distance
  )

  # Return results
  list(
    curves = matrix(unlist(curves), nrow = 100, byrow = TRUE),
    noisy_curves = matrix(unlist(noisy_curves), nrow = 100, byrow = TRUE),
    noisy_curves_miss = noisy_curves_miss,
    response = response,
    grid = t,
    beta = beta,
    stochastic_components = stochastic_components
  )
}

# # Example usage
# set.seed(123)
#
# # Generate data with random component
# data_random <- generate_1d_po_functional_data(n = 5, beta_type = "gaussian")
# print("Random stochastic components:")
# print(data_random$stochastic_components)
#
# # Visualize example curves
# par(mfrow = c(1, 3))
# plot(data_random$grid, data_random$curves[[2]], type = "l",
#      main = "True Curve", xlab = "t", ylab = "y")
# plot(data_random$grid, data_random$noisy_curves[[2]], type = "l",
#      main = "Noisy Curve", xlab = "t", ylab = "y")
# plot(data_random$grid, data_random$beta, type = "l",
#      main = "Beta Coefficient", xlab = "t", ylab = "y")
#
# plot(data_random$grid, data_random$noisy_curves_miss[[1]][[1]], type = "l",
#      main = "Miss Curve", xlab = "t", ylab = "y")
