#' Data generator function for the variable domain case
#'
#' Generates a variable domain functional regression model
#'
#' @param N Number of subjects.
#' @param J Number of maximum observations per subject.
#' @param nsims Number of simulations per the simulation study.
#' @param aligned If the data that will be generated is aligned or not.
#' @param multivariate If TRUE, the data is generated with 2 functional variables.
#' @param beta_index Index for the beta.
#' @param Rsq Variance of the model.
#' @param use_x If the data is generated with x.
#' @param use_f If the data is generated with f.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item y: \code{vector} of length N containing the response variable.
#'   \item X_s: \code{matrix} of non-noisy functional data for the first functional covariate.
#'   \item X_se: \code{matrix} of noisy functional data for the first functional covariate
#'   \item Y_s: \code{matrix} of non-noisy functional data for the second functional covariate (if multivariate).
#'   \item Y_se: \code{matrix} of noisy functional data for the second covariate (if multivariate).
#'   \item x1: \code{vector} of length N containing the non-functional covariate (if use_x is TRUE).
#'   \item x2: \code{vector} of length N containing the observed values of the smooth term (if use_f is TRUE).
#'   \item smooth_term: \code{vector} of length N containing a smooth term (if use_f is TRUE).
#'   \item Beta: \code{array} containing the true functional coefficients.
#' }
#'
#' @examples
#' # Basic usage with default parameters
#' sim_data <- data_generator_vd()
#'
#' # Generate data with non-aligned domains
#' non_aligned_data <- data_generator_vd(N = 150, J = 120, aligned = FALSE)
#'
#' # Generate multivariate functional data
#' multivariate_data <- data_generator_vd(N = 200, J = 100, multivariate = TRUE)
#'
#' # Generate data with non-functional covariates and smooth term
#' complex_data <- data_generator_vd(
#'   N = 100,
#'   J = 150,
#'   use_x = TRUE,
#'   use_f = TRUE
#' )
#'
#' # Generate data with a different beta function and R-squared value
#' custom_beta_data <- data_generator_vd(
#'   N = 80,
#'   J = 80,
#'   beta_index = 2,
#'   Rsq = 0.8
#' )
#'
#' # Access components of the generated data
#' y <- sim_data$y # Response variable
#' X_s <- sim_data$X_s # Noise-free functional covariate
#' X_se <- sim_data$X_se # Noisy functional covariate
#'
#' @export
data_generator_vd <- function(
    N = 100,
    J = 100,
    nsims = 1,
    Rsq = 0.95,
    aligned = TRUE,
    multivariate = FALSE,
    beta_index = 1,
    use_x = FALSE,
    use_f = FALSE
) {
  if (!(beta_index %in% c(1, 2))) {
    stop("'beta_index' could only be 1 or 2", call. = FALSE)
  }

  for (iter in 1:nsims) {

    if (aligned) {
      # Generating the domain for all subject with a minimum of 10 observations (min = 10)
      M <- round(stats::runif(N, min = 10, max = J), digits = 0)

      M <- sort(M) # We can sort the data without loss of generality
    } else {
      M <- cbind(
        round(stats::runif(N, min = 1, max = (J / 2) - 5), digits = 0),
        round(stats::runif(N, min = (J / 2) + 5, max = J), digits = 0)
      )

      M_diff <- M[, 2] - M[, 1] + 1
    }

    maxM <- max(M)
    t <- 1:maxM

    # Here we generate the functional data
    X_s <- matrix(NA, N, maxM) # NOT NOISY
    X_se <- matrix(NA, N, maxM) # NOISY
    Y_s <- matrix(NA, N, maxM) # NOT NOISY
    Y_se <- matrix(NA, N, maxM) # NOISY

    for (i in 1:N) {
      u1 <- stats::rnorm(1)

      temp <- matrix(NA, 10, maxM)

      for (k in 1:10) {
        v_i1 <- stats::rnorm(1, 0, 4 / k^2)
        v_i2 <- stats::rnorm(1, 0, 4 / k^2)

        if (aligned) {
          temp[k, 1:M[i]] <-
            v_i1 * sin(2 * pi * k * (1:M[i]) / J) + v_i2 * cos(2 * pi * k * (1:M[i]) / J)
        } else {
          temp[k, (M[i, 1]:M[i, 2])] <-
            v_i1 * sin(2 * pi * k * (M[i, 1]:M[i, 2]) / J) + v_i2 * cos(2 * pi * k * (M[i, 1]:M[i, 2]) / J)
        }
      }

      B <- apply(temp, 2, sum)


      B <- B + u1
      B2 <- B + stats::rnorm(1, sd = 0.02) + (t / 10)

      aux <- stats::var(B, na.rm = TRUE)

      X_s[i, ] <- B
      X_se[i, ] <- B + stats::rnorm(maxM, 0, sqrt(aux / 8)) # WE ADD NOISE
      Y_s[i, ] <- B2
      Y_se[i, ] <- B2 + stats::rnorm(maxM, 0, sqrt(aux / 8)) # WE ADD NOISE
    }

    Beta <- array(dim = c(N, maxM, 4))
    nu <- rep(0, N)
    y <- rep(0, N)
    x1 <- stats::rnorm(N)
    x2 <- stats::runif(N)
    f1 <- function(x) 2 * sin(pi * x)
    f2 <- function(x) 3.5 * cos(pi * x)

    for (i in 1:N) {
      # Computing the true functional coefficients
      if (aligned) {
        Beta[i, 1:(M[i]), 1] <- ((10 * t[1:(M[i])] / M[i]) - 5) / 10
        Beta[i, 1:(M[i]), 2] <- ((1 - (2 * M[i] / maxM)) * (5 - 40 * ((t[1:(M[i])] / M[i]) - 0.5)^2)) / 10

        if (multivariate) {
          nu[i] <- sum(X_s[i, ] * Beta[i, , beta_index], na.rm = TRUE) / (M[i]) + sum(Y_s[i, ] * Beta[i, , 2], na.rm = TRUE) / (M[i])
        } else {
          nu[i] <- sum(X_s[i, ] * Beta[i, , beta_index], na.rm = TRUE) / (M[i]) # NOT NOISY
        }
      } else {
        Beta[i, (M[i, 1]:M[i, 2]), 1] <- ((10 * t[(M[i, 1]:M[i, 2])] / M_diff[i]) - 5) / 10
        Beta[i, (M[i, 1]:M[i, 2]), 2] <- ((1 - (2 * M_diff[i] / maxM)) * (5 - 40 * ((t[(M[i, 1]:M[i, 2])] / M_diff[i]) - 0.5)^2)) / 10

        if (multivariate) {
          nu[i] <- sum(X_s[i, ] * Beta[i, , beta_index], na.rm = TRUE) / (M_diff[i]) + sum(Y_s[i, ] * Beta[i, , 2], na.rm = TRUE) / (M_diff[i])
        } else {
          nu[i] <- sum(X_s[i, ] * Beta[i, , beta_index], na.rm = TRUE) / (M_diff[i]) # NOT NOISY
        }
      }
    }

    smooth_term <- f1(x2)
    nu <- if (use_f) nu + smooth_term else nu
    var_e <- (1 / Rsq - 1) * stats::var(nu)
    y <- nu + stats::rnorm(N, sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL
    y <- if (use_x) y + x1 else y
  }


  data <- list(y = y)

  data[["X_s"]] <- X_s
  data[["X_se"]] <- X_se
  data[["Y_s"]] <- Y_s
  data[["Y_se"]] <- Y_se
  data[["x1"]] <- x1
  data[["x2"]] <- x2
  data[["smooth_term"]] <- smooth_term
  data[["Beta"]] <- Beta

  data
}


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
      curve_noisy = true_curve + stats::rnorm(length(t), 0, noise_sd)
    )
  }


  # Generate coefficient function
  beta <- if (beta_type == "sin") {
    sin(2 * pi * t) * cos(pi * t)
  } else {
    stats::dnorm(t, mean = 0.5, sd = 0.15) / stats::dnorm(0.5, mean = 0.5, sd = 0.15)
  }

  # Generate data and compute response
  for (i in 1:n) {
    # Use fixed value if provided, otherwise generate random component
    stochastic_components[[i]] <- stats::rnorm(3, 0, 0.2)

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
  var_e <- (1/rsq - 1) * stats::var(nu)
  response <- nu + stats::rnorm(n, 0, sqrt(var_e))

  # Add missing values
  noisy_curves_miss <- add_miss1d(
    noisy_curves,
    n_missing = n_missing,
    min_distance = min_distance
  )

  miss_points <- noisy_curves_miss[["miss_points"]]
  missing_points <- noisy_curves_miss[["missing_points"]]
  X_miss <- matrix(
    unlist(noisy_curves_miss[["X_miss"]]),
    nrow = length(noisy_curves_miss[["X_miss"]]),
    ncol = length(noisy_curves_miss[["X_miss"]][[1]]),
    byrow = TRUE
  )

  # Return results
  list(
    curves = matrix(unlist(curves), nrow = 100, byrow = TRUE),
    noisy_curves = matrix(unlist(noisy_curves), nrow = 100, byrow = TRUE),
    noisy_curves_miss = X_miss,
    miss_points = miss_points,
    missing_points = missing_points,
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
#' @param n_missing Number of holes in every curve
#' @param min_distance_x Length of the holes in the x axis
#' @param min_distance_y Length of the holes in the y axis
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
    n = 20,
    grid_x = 20,
    grid_y = 20,
    noise_sd = 0.015,
    rsq = 0.95,
    beta_type = c("saddle", "exp"),
    response_type = c("gaussian","binomial"),
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
  stochastic_components <- matrix(
    nrow = n,
    ncol = 2,
    dimnames = list(NULL, c("a1", "a2"))
  )
  nu <- numeric(n)

  # Generate data and compute response
  for (i in 1:n) {
    # Use fixed values if provided, otherwise generate random components
    if (!is.null(a1)) {
      stochastic_components[i, ] <- c(a1, a2)
    } else {
      stochastic_components[i, ] <- c(
        stats::rnorm(1, 0, 0.2),  # a1
        stats::rnorm(1, 0, 0.2)   # a2
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

  response <- if (response_type == "gaussian") {
    # Generate response with desired R-squared
    var_e <- (1/rsq - 1) * stats::var(nu)
    response <- nu + stats::rnorm(n, 0, sqrt(var_e))

  } else {
    # BINOMIAL MODEL
    response <- rbinom(N,1,(exp(nu)/(1+exp(nu))))

  }



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
    noisy_surfaces_miss = noisy_surfaces_miss[[1]],
    miss_points = noisy_surfaces_miss[[2]],
    missing_points = noisy_surfaces_miss[[3]],
    response = response,
    grid_x = x,
    grid_y = y,
    beta = beta_surface,
    stochastic_components = stochastic_components
  )
}

#' Generate saddle-shaped coefficient surface
#'
#' @noRd
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
#' @noRd
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

#' Helper function to generate surface with given parameters
#'
#' @noRd
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
    DATA_N = true_surface + matrix(
      stats::rnorm(length(x) * length(y), 0, noise_sd),
      length(x),
      length(y)
    )
  )
}

#' Setup Simpson's integration weights
#'
#' @noRd
setup_simpson_weights <- function(n_x, n_y, h_x, h_y) {
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
