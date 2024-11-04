#' Set missing values to the surfaces
#'
#' @param X List of surfaces.
#' @param n_missing Number of holes in every surface.
#' @param min_distance_x Length of the holes in the x-axis.
#' @param min_distance_y Length of the holes in the y-axis.
#'
#' @return List of surfaces with missing values. The difference between
#' miss_points and missing_points is the format in which the data is
#' presented.
#'
#' @noRd
add_miss2 <- function(X, n_missing = 1, min_distance_x = 9, min_distance_y = 9) {
  N <- length(X)
  missing_points <- miss_points <- vector(mode = "list", length = N)

  x_b <- nrow(X[[1]])
  y_b <- ncol(X[[1]])

  for (i in 2:N) {
    if (nrow(X[[i]]) != x_b || ncol(X[[i]]) != y_b) {
      stop("The dimension of all the surfaces must be the same.", call. = FALSE)
    }
  }

  for (j in seq_along(missing_points)) {

    if (n_missing>=1) {

      x_missing <- sort(sample(1:(x_b - min_distance_x), n_missing))
      y_missing <- sort(sample(1:(y_b - min_distance_y), n_missing))

      x_pos <- x_missing:(x_missing + min_distance_x - 1)
      y_pos <- y_missing:(y_missing + min_distance_y - 1)

      missing_points[[j]] <- expand.grid(x_pos, y_pos)

      for (add_miss in 1:nrow(missing_points[[j]])) {
        X[[j]][missing_points[[j]][add_miss, 1], missing_points[[j]][add_miss, 2]] <-
          NA
      }
    }


    miss_points[[j]] <- vector(mode = "list", length = ncol(X[[j]]))

    for (i in 1:ncol(X[[j]])) {
      miss_spots <- NULL

      for (j_row in 1:nrow(X[[j]])) {
        if (is.na(X[[j]][j_row, i])) {
          miss_spots <- c(miss_spots, j_row)
        }
      }

      if (!is.null(miss_spots)) {
        miss_points[[j]][[i]] <- miss_spots
      }
    }
  }

  list(
    X_miss = X,
    miss_points = miss_points,
    missing_points = missing_points
  )
}

#' Generate 2D functional data for simulation studies
#'
#' Creates synthetic 2D functional data with optional noise components and different
#' coefficient patterns. Uses Simpson's rule for accurate integration.
#'
#' @param n Number of samples to generate
#' @param grid_x Number of points in x-axis grid. Default is 20
#' @param grid_y Number of points in y-axis grid. Default is 20
#' @param noise_sd Standard deviation of measurement noise. Default is 0.015
#' @param rsq Desired R-squared value for the response. Default is 0.95
#' @param beta_type Type of coefficient surface ("saddle" or "exp"). Default is "saddle"
#' @param a1 Optional fixed value for first stochastic component. If provided, a2 must also be provided
#' @param a2 Optional fixed value for second stochastic component. If provided, a1 must also be provided
#' @param sub_response Number of intervals for Simpson integration. Default is 50
#'
#' @return A list containing:
#' \itemize{
#'   \item surfaces: List of n true (noiseless) surfaces
#'   \item noisy_surfaces: List of n observed (noisy) surfaces
#'   \item response: Vector of n response values
#'   \item grid_x: x-axis grid points
#'   \item grid_y: y-axis grid points
#'   \item beta: True coefficient surface
#'   \item stochastic_components: Matrix of a1 and a2 values used for each surface
#' }
#'
generate_2d_po_functional_data <- function(
  n = 100,
  grid_x = 20,
  grid_y = 20,
  noise_sd = 0.015,
  rsq = 0.95,
  beta_type = c("saddle", "exp"),
  a1 = NULL,
  a2 = NULL,
  sub_response = 50,
  n_missing = 1,
  min_distance_x = NULL,
  min_distance_y = NULL
) {

  beta_type <- match.arg(beta_type)

  # Validate a1 and a2 parameters
  if (!is.null(a1) && is.null(a2) || is.null(a1) && !is.null(a2)) {
    stop("Both a1 and a2 must be provided if one is specified", call. = FALSE)
  }

  if (!is.null(a1)) {
    if (length(a1) != 1 || length(a2) != 1 || !is.numeric(a1) || !is.numeric(a2)) {
      stop("a1 and a2 must be single numeric values", call. = FALSE)
    }
  }

  if (!is.null(min_distance_x) && is.null(min_distance_y) || is.null(min_distance_x) && !is.null(min_distance_y)) {
    stop("Both min_distance_x and min_distance_y must be provided if one is specified", call. = FALSE)
  }

  if (is.null(min_distance_x)) {
    min_distance_x <- round(1/5*grid_x)
  }

  if (is.null(min_distance_y)) {
    min_distance_y <- round(1/5*grid_y)
  }

  # Create grid points
  x <- seq(0, 1, length.out = grid_x)
  y <- seq(0, 1, length.out = grid_y)

  # Initialize storage
  surfaces <- vector("list", n)
  noisy_surfaces <- vector("list", n)
  stochastic_components <- matrix(nrow = n, ncol = 2,
                                  dimnames = list(NULL, c("a1", "a2")))
  nu <- numeric(n)

  # Helper function to generate surface with given parameters
  generate_surface <- function(x, y, a1, a2, noise_sd) {
    true_surface <- matrix(nrow = length(x), ncol = length(y))
    for (i in seq_along(x)) {
      for (j in seq_along(y)) {
        true_surface[i, j] <- a1 * cos(2 * pi * x[i]) +
          a2 * cos(2 * pi * y[j]) + 1
      }
    }
    list(
      DATA_T = true_surface,
      DATA_N = true_surface + matrix(rnorm(length(x) * length(y), 0, noise_sd),
                                     length(x), length(y))
    )
  }

  # Setup Simpson's integration weights
  setup_simpson_weights <- function(n_x, n_y, h_x, h_y) {
    # Initialize weight vector
    W_delta <- array(dim = (n_x + 1) * (n_y + 1))

    # Create Simpson weights for x direction
    simp_w_x <- rep(1, n_x + 1)
    simp_w_x[seq(2, n_x - 1, 2)] <- 4
    simp_w_x[seq(3, n_x - 1, 2)] <- 2

    # Combine with h_x/3
    W_x <- (h_x/3) * simp_w_x

    # Create full weight matrix
    for (j in 1:(n_y + 1)) {
      start_idx <- ((n_x + 1) * (j - 1) + 1)
      end_idx <- ((n_x + 1) * j)

      if (j == 1 || j == (n_y + 1)) {
        W_delta[start_idx:end_idx] <- (h_y/3) * W_x
      } else if (j %% 2 == 0) {
        W_delta[start_idx:end_idx] <- (4 * h_y/3) * W_x
      } else {
        W_delta[start_idx:end_idx] <- (2 * h_y/3) * W_x
      }
    }

    W_delta
  }

  # Generate data and compute response
  for (i in 1:n) {
    # Use fixed values if provided, otherwise generate random components
    if (!is.null(a1)) {
      stochastic_components[i, ] <- c(a1, a2)
    } else {
      stochastic_components[i, ] <- c(
        rnorm(1, 0, 0.2),  # a1
        rnorm(1, 0, 0.2)   # a2
      )
    }

    # Generate surface
    surface_data <- generate_surface(x, y,
                                     stochastic_components[i, 1],
                                     stochastic_components[i, 2],
                                     noise_sd)

    surfaces[[i]] <- surface_data$DATA_T
    noisy_surfaces[[i]] <- surface_data$DATA_N

    # Prepare integration grid
    n_x <- m_y <- 2 * sub_response

    x_fine <- seq(min(x), max(x), length.out = n_x + 1)
    y_fine <- seq(min(y), max(y), length.out = m_y + 1)

    h_x <- (max(x_fine) - min(x_fine)) / n_x
    h_y <- (max(y_fine) - min(y_fine)) / m_y

    # Generate finer surface for integration
    surface_fine <- generate_surface(x_fine, y_fine,
                                     stochastic_components[i, 1],
                                     stochastic_components[i, 2],
                                     0)$DATA_T

    # Generate beta surface on fine grid
    beta_fine <- if (beta_type == "saddle") {
      generate_saddle_surface(x_fine, y_fine)
    } else {
      generate_exp_surface(x_fine, y_fine)
    }

    # Get integration weights
    W_delta <- setup_simpson_weights(n_x, m_y, h_x, h_y)

    # Compute double integral using Simpson's rule
    nu[i] <- as.double(t(as.vector(surface_fine)) %*%
                         diag(W_delta) %*%
                         as.vector(beta_fine))
  }

  # Generate response with desired R-squared
  var_e <- (1/rsq - 1) * var(nu)
  response <- nu + rnorm(n, 0, sqrt(var_e))

  # Generate beta on original grid for output
  beta_surface <- if (beta_type == "saddle") {
    generate_saddle_surface(x, y)
  } else {
    generate_exp_surface(x, y)
  }

  noisy_surfaces_miss <- add_miss2(
    noisy_surfaces,
    n_missing,
    min_distance_x,
    min_distance_y
  )

  # Return results
  list(
    surfaces = surfaces,
    noisy_surfaces = noisy_surfaces,
    noisy_surfaces_miss = noisy_surfaces_miss,
    response = response,
    grid_x = x,
    grid_y = y,
    beta = beta_surface,
    stochastic_components = stochastic_components
  )
}

#' Generate saddle-shaped coefficient surface
#'
#' @param x Grid points for x-axis
#' @param y Grid points for y-axis
#' @return Matrix containing the coefficient surface values
#'
generate_saddle_surface <- function(x, y) {
  surface <- matrix(nrow = length(x), ncol = length(y))

  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      surface[i, j] <- ((4 * x[i] - 2)^3 -
                          3 * (4 * x[i] - 2) * (4 * y[j] - 2)^2) / 10
    }
  }
  surface
}

#' Generate exponential coefficient surface
#'
#' @param x Grid points for x-axis
#' @param y Grid points for y-axis
#' @return Matrix containing the coefficient surface values
#'
generate_exp_surface <- function(x, y) {
  surface <- matrix(nrow = length(x), ncol = length(y))

  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      surface[i, j] <- (5 * exp(-8 * ((x[i] - 0.75)^2 + (y[j] - 0.75)^2)) +
                          5 * exp(-8 * ((x[i] - 0.1)^2 + (y[j] - 0.1)^2))) / 10
    }
  }
  surface
}

# # Example usage
# set.seed(123)
#
# # Generate data with random a1 and a2
# data_random <- generate_2d_po_functional_data(n = 5)
# print("Random stochastic components:")
# print(data_random$stochastic_components)
#
# # Generate data with fixed a1 and a2
# data_fixed <- generate_2d_po_functional_data(n = 5, a1 = 0.5, a2 = -0.3)
# print("\nFixed stochastic components:")
# print(data_fixed$stochastic_components)
#
# data <- generate_2d_po_functional_data()
# plotly::plot_ly(z = data$noisy_surfaces_miss$X_miss[[1]], type = "surface")
#
# # Visualize example surfaces
# par(mfrow = c(1, 3))
# image(data_fixed$grid_x, data_fixed$grid_y, data_fixed$surfaces[[1]],
#       main = "True Surface (Fixed a1,a2)", xlab = "x", ylab = "y")
# image(data_fixed$grid_x, data_fixed$grid_y, data_fixed$noisy_surfaces[[1]],
#       main = "Noisy Surface", xlab = "x", ylab = "y")
# image(data_fixed$grid_x, data_fixed$grid_y, data_fixed$beta,
#       main = "Beta Coefficient", xlab = "x", ylab = "y")
