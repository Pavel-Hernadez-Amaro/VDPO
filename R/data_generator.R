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
    use_f = FALSE) {
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

#' Add missing values to curves
#'
#' @param X List of curves
#' @param n_missing Number of holes in every curve
#' @param min_distance Length of the holes
#'
#' @return List containing curves with missing values and missing points information
#'
#' @noRd
add_miss1d_end <- function(X, n_missing = 1, min_distance = 5) {
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
      if (j<floor(length(seq_along(missing_points))/2)) {

        x_missing <- rep(1, n_missing)

      }else{
        x_missing <- rep(n_points - min_distance + 1, n_missing)
        }


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
#' Generate 1D Functional Data for Simulation Studies
#'
#' Creates synthetic 1D functional data with optional noise components and different
#' coefficient patterns. Uses the trapezoidal rule for numerical integration.
#'
#' @param n Number of samples to generate. Default is 100.
#' @param grid_points Number of points in the grid. Default is 100.
#' @param noise_sd Standard deviation of measurement noise. Default is 0.015.
#' @param rsq Desired R-squared value for the response. Default is 0.95.
#' @param beta_type Type of coefficient function ("sin" or "gaussian"). Default is "sin".
#' @param n_missing Number of missing segments per curve. Default is 1.
#' @param min_distance Minimum length of missing segments. Default is NULL (auto-calculated).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{curves}: Matrix of \code{n} true (noiseless) curves, each as a row.
#'   \item \code{noisy_curves}: Matrix of \code{n} observed (noisy) curves, each as a row.
#'   \item \code{noisy_curves_miss}: Matrix of noisy curves with missing values.
#'   \item \code{miss_points}: Indices of the missing segments in the noisy curves.
#'   \item \code{missing_points}: Details of the missing segments for each curve.
#'   \item \code{response}: Vector of \code{n} response values.
#'   \item \code{grid}: Grid points on which the curves are defined.
#'   \item \code{beta}: Coefficient function applied to the curves.
#'   \item \code{stochastic_components}: List of stochastic coefficients used for each curve.
#' }
#'
#' @examples
#' # Generate basic 1D functional data with default parameters
#' data <- data_generator_po_1d(n = 10)
#'
#' # Generate data with a Gaussian-shaped coefficient function
#' data <- data_generator_po_1d(n = 2, beta_type = "gaussian")
#'
#' # Generate data with higher grid resolution
#' data <- data_generator_po_1d(n = 2, grid_points = 200)
#'
#' # Generate data with larger measurement noise
#' data <- data_generator_po_1d(n = 2, noise_sd = 0.05)
#'
#' # Introduce missing segments in the curves
#' data <- data_generator_po_1d(n = 2, n_missing = 3, min_distance = 10)
#'
#' # Generate data with low R-squared value
#' data <- data_generator_po_1d(n = 2, rsq = 0.8)
#'
#' @export
data_generator_po_1d_old <- function(
    n = 100,
    grid_points = 100,
    noise_sd = 0.25,
    center = TRUE,
    rsq = 0.95,
    mu = 0.1,
    beta_type = c("sin", "trig", "exp", "linear","quadratic", "cubic", "Wang"),
    beta_type_2 = c("sin", "trig", "exp", "linear","quadratic", "cubic", "Wang"),
    univariate = TRUE,
    response_type = c("gaussian","binomial"),
    linear_predictor = c("integral","linear"),
    n_missing = 1,
    min_distance = NULL) {

  beta_type <- match.arg(beta_type)
  beta_type_2 <- match.arg(beta_type_2)

  if (is.null(min_distance)) {
    min_distance <- round(1 / 5 * grid_points)
  }

  # Create grid points
  t <- seq(0, 1, length.out = grid_points)

  # Initialize storage
  curves <- vector("list", n)
  noisy_curves <- vector("list", n)
  stochastic_components <- vector("list", n)
  if (!univariate) {
    curves_2 <- vector("list", n)
    noisy_curves_2 <- vector("list", n)
    stochastic_components_2 <- vector("list", n)
  }
  nu <- numeric(n)

  # # Helper function to generate curve with given parameters
  # generate_curve_1 <- function(t, a1, a2, b1, noise_sd) {
  #   true_curve <- a1 * sin(2 * pi * t) +
  #     a2 * cos(4 * pi * t) +
  #     b1 * t
  #
  #   if (center == TRUE) {
  #     sol=list(
  #       curve_true = true_curve - mean(true_curve),
  #       curve_noisy = true_curve - mean(true_curve) + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
  #     )
  #   }else{
  #     sol=list(
  #       curve_true = true_curve,
  #       curve_noisy = true_curve + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
  #     )
  #   }
  #   sol
  #   }
  #
  # generate_curve_2 <- function(t, a1, a2, b1, noise_sd) {
  #   true_curve <- a1 * sin(2 * pi * t) +
  #     a2 * cos(4 * pi * t) +
  #     exp(-t) +
  #     1
  #
  #   if (center == TRUE) {
  #     sol=list(
  #       curve_true = true_curve - mean(true_curve),
  #       curve_noisy = true_curve - mean(true_curve) + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
  #     )
  #   }else{
  #     sol=list(
  #       curve_true = true_curve,
  #       curve_noisy = true_curve + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
  #     )
  #   }
  #   sol
  #   }

  # Helper function to generate curve with given parameters
  generate_curve_1 <- function(t, a1, a2, a3, noise_sd) {
    true_curve <- a1 * cos(1 * pi * t) +
                  a2 * cos(3 * pi * t) +
                  a3 * cos(5 * pi * t)


    if (center == TRUE) {
      sol=list(
        curve_true = true_curve - mean(true_curve),
        curve_noisy = true_curve - mean(true_curve) + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
      )
    }else{
      sol=list(
        curve_true = true_curve,
        curve_noisy = true_curve + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
      )
    }
    sol
  }

  generate_curve_2 <- function(t, a1, a2, a3, noise_sd) {
    true_curve <- a1 * sin(1 * pi * t) +
                  a2 * sin(3 * pi * t) +
                  a3 * sin(5 * pi * t)


    if (center == TRUE) {
      sol=list(
        curve_true = true_curve - mean(true_curve),
        curve_noisy = true_curve - mean(true_curve) + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
      )
    }else{
      sol=list(
        curve_true = true_curve,
        curve_noisy = true_curve + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
      )
    }
    sol
  }

  # Function to generate beta coefficients as in your original code
  generate_beta <- function(beta_type, t) {
    if (beta_type == "sin") {
      2*sin(0.5*pi*(t))+4*sin(1.5*pi*(t))+5*sin(2.5*pi*(t))
    } else if(beta_type == "trig") {
      0.5 * sin(2 * pi * t) + 2 * cos(pi * t)
    } else if(beta_type == "cubic") {
      0.5*t^3
    } else if(beta_type == "exp") {
      0.5 * exp(t)
    } else if(beta_type == "linear") {
      2*t
    } else if(beta_type == "Wang") {
      1+3*sqrt(2)*cos(2 * pi * t)
    } else { # quadratic
      -t^2
    }
  }

  # Function to scale a vector to a desired range
  scale_to_range <- function(vec) {

    # Check if beta2 can be negative
    if(min(vec) < 0){
      new_min = -1
    }else{
      new_min = 0
      }

    new_max = 1
    old_min <- min(vec)
    old_max <- max(vec)



    # Check if vector is constant
    if (old_min == old_max) {
      return(rep((new_min + new_max) / 2, length(vec)))
    }

    # Scale
    scaled <- (vec - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
    return(scaled)
  }

  # Create grid
  t <- seq(0, 1, length.out = 100)

  # Generate beta coefficients with your original functions
  beta <- generate_beta(beta_type, t)
  beta_scaled <- scale_to_range(beta)

  # # Generate coefficient function
  # beta <- if (beta_type == "sin") {
  #   2*sin(0.5*pi*(t))+4*sin(1.5*pi*(t))+5*sin(2.5*pi*(t))
  # } else if(beta_type == "trig") {
  #   0.5 * sin(2 * pi * t) + 2 * cos(pi * t)
  #   }else if(beta_type == "cubic") {
  #     0.5*t^3
  #   }else if(beta_type == "exp") {
  #     0.5 * exp(t)
  # }else if(beta_type == "linear"){
  #   2*t
  # }else if(beta_type == "Wang"){
  #   1+3*sqrt(2)*cos(2 * pi * t)
  # }else{
  #   -t^2
  # }

  # Generate data and compute response
  for (i in 1:n) {
    # Use fixed value if provided, otherwise generate random component
    stochastic_components[[i]] <- stats::rnorm(3, 0, 0.2)

    # Generate curve
    curve_data <- generate_curve_1(
      t,
      stochastic_components[[i]][1],
      stochastic_components[[i]][2],
      stochastic_components[[i]][3],
      noise_sd
    )

    curves[[i]] <- curve_data$curve_true
    noisy_curves[[i]] <- curve_data$curve_noisy

    if (univariate) {

      if (linear_predictor=="integral") {
      integrand <- curves[[i]] * beta_scaled
      nu[i] <- mu + (0.5 * sum(diff(t) * (integrand[-1] + integrand[-length(integrand)])))
      }else{
        nu[i] <- mu + ((curves[[i]] %*% beta_scaled) / grid_points)
      }

    }else{

      # Use fixed value if provided, otherwise generate random component
      stochastic_components_2[[i]] <- stats::rnorm(3, 0, 0.2)

      # Generate curve
      curve_data_2 <- generate_curve_2(
        t,
        stochastic_components_2[[i]][1],
        stochastic_components_2[[i]][2],
        stochastic_components_2[[i]][3],
        noise_sd
      )

      curves_2[[i]] <- curve_data_2$curve_true
      noisy_curves_2[[i]] <- curve_data_2$curve_noisy

      # Generate coefficient function
      # CHANGED BETA TYPES SO IT COULD HAVE DIFFERENT ONES FOR EACH VARIABLE
      # beta_2 <- if (beta_type_2 == "sin") {
      #   2*sin(0.5*pi*(t))+4*sin(1.5*pi*(t))+5*sin(2.5*pi*(t))
      # }else if(beta_type_2 == "trig") {
      #   0.5 * sin(2 * pi * t) + 2 * cos(pi * t)
      # }else if(beta_type_2 == "cubic") {
      #   -0.5*t^3
      # }else if(beta_type_2 == "exp") {
      #   0.5 * exp(t)
      # }else if(beta_type_2 == "linear"){
      #   -2*t
      # }else if(beta_type_2 == "Wang"){
      #  1+3*sqrt(2)*cos(2 * pi * t)
      # }else{
      #   t^2
      # }

      beta_2 <- generate_beta(beta_type_2, t)
      beta_2_scaled <- scale_to_range(beta_2)


      integrand <- curves[[i]] * beta_scaled
      integrand_2 <- curves_2[[i]] * beta_2_scaled

      if (linear_predictor=="integral") {
        nu[i] <- mu + ( 0.5 * sum(diff(t) * (integrand[-1] + integrand[-length(integrand)])) +
          0.5 * sum(diff(t) * (integrand_2[-1] + integrand_2[-length(integrand_2)])))
      }else{
        nu[i] <- mu + (((curves[[i]] %*% beta_scaled) + (curves_2[[i]] %*% beta_2_scaled)) / grid_points)
      }

    }
    # Compute integral using trapezoidal rule
  }

  if (response_type=="gaussian") {
    # Generate response with desired R-squared
    var_e <- (1 / rsq - 1) * stats::var(nu)
    print(var_e)
    response <- nu + stats::rnorm(n, 0, sqrt(var_e))
  }else{
    response <- stats::rbinom(n, 1, (exp(nu) / (1 + exp(nu))))
    # if (sum(response)>=65 || sum(response)<=35) {
    #   zeros <- round(n/2)
    #   response <- c(rep(0,zeros),rep(1,n-zeros))
    # }
  }

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

  if (!univariate) {
    noisy_curves_miss_2 <- add_miss1d(
      noisy_curves_2,
      n_missing = n_missing,
      min_distance = min_distance
    )

    miss_points_2 <- noisy_curves_miss_2[["miss_points"]]
    missing_points_2 <- noisy_curves_miss_2[["missing_points"]]
    X_miss_2 <- matrix(
      unlist(noisy_curves_miss_2[["X_miss"]]),
      nrow = length(noisy_curves_miss_2[["X_miss"]]),
      ncol = length(noisy_curves_miss_2[["X_miss"]][[1]]),
      byrow = TRUE
    )
    }



  # Return results
  return(if (univariate) {
    list(
      curves = matrix(unlist(curves), nrow = n, byrow = TRUE),
      noisy_curves = matrix(unlist(noisy_curves), nrow = n, byrow = TRUE),
      noisy_curves_miss = X_miss,
      miss_points = miss_points,
      missing_points = missing_points,
      response = response,
      grid = t,
      beta = beta,
      beta_scaled = beta_scaled,
      stochastic_components = stochastic_components
    )

  }else{
    list(
      curves_1 = matrix(unlist(curves), nrow = n, byrow = TRUE),
      noisy_curves_1 = matrix(unlist(noisy_curves), nrow = n, byrow = TRUE),
      noisy_curves_miss_1 = X_miss,
      miss_points_1 = miss_points,
      missing_points_1 = missing_points,
      curves_2 = matrix(unlist(curves_2), nrow = n, byrow = TRUE),
      noisy_curves_2 = matrix(unlist(noisy_curves_2), nrow = n, byrow = TRUE),
      noisy_curves_miss_2 = X_miss_2,
      miss_points_2 = miss_points_2,
      missing_points_2 = missing_points_2,
      beta = beta,
      beta_2 = beta_2,
      beta_scaled = beta_scaled,
      beta_2_scaled = beta_2_scaled,
      stochastic_components = stochastic_components,
      stochastic_components_2 = stochastic_components_2,
      grid = t,
      response = response
    )

  }
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
#' Generate 1D Functional Data for Simulation Studies
#'
#' Creates synthetic 1D functional data with optional noise components and different
#' coefficient patterns. Uses the trapezoidal rule for numerical integration.
#'
#' @param n Number of samples to generate. Default is 100.
#' @param grid_points Number of points in the grid. Default is 100.
#' @param noise_sd Standard deviation of measurement noise. Default is 0.015.
#' @param rsq Desired R-squared value for the response. Default is 0.95.
#' @param beta_type Type of coefficient function ("sin" or "gaussian"). Default is "sin".
#' @param n_missing Number of missing segments per curve. Default is 1.
#' @param min_distance Minimum length of missing segments. Default is NULL (auto-calculated).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{curves}: Matrix of \code{n} true (noiseless) curves, each as a row.
#'   \item \code{noisy_curves}: Matrix of \code{n} observed (noisy) curves, each as a row.
#'   \item \code{noisy_curves_miss}: Matrix of noisy curves with missing values.
#'   \item \code{miss_points}: Indices of the missing segments in the noisy curves.
#'   \item \code{missing_points}: Details of the missing segments for each curve.
#'   \item \code{response}: Vector of \code{n} response values.
#'   \item \code{grid}: Grid points on which the curves are defined.
#'   \item \code{beta}: Coefficient function applied to the curves.
#'   \item \code{stochastic_components}: List of stochastic coefficients used for each curve.
#' }
#'
#' @examples
#' # Generate basic 1D functional data with default parameters
#' data <- data_generator_po_1d(n = 10)
#'
#' # Generate data with a Gaussian-shaped coefficient function
#' data <- data_generator_po_1d(n = 2, beta_type = "gaussian")
#'
#' # Generate data with higher grid resolution
#' data <- data_generator_po_1d(n = 2, grid_points = 200)
#'
#' # Generate data with larger measurement noise
#' data <- data_generator_po_1d(n = 2, noise_sd = 0.05)
#'
#' # Introduce missing segments in the curves
#' data <- data_generator_po_1d(n = 2, n_missing = 3, min_distance = 10)
#'
#' # Generate data with low R-squared value
#' data <- data_generator_po_1d(n = 2, rsq = 0.8)
#'
#' @export
data_generator_po_1d <- function(
    n = 100,
    grid_points = 100,
    noise_sd = 0.25,
    center = TRUE,
    rsq = 0.95,
    mu = 0.1,
    univariate = TRUE,
    response_type = c("gaussian","binomial"),
    linear_predictor = c("rectangular", "trapezoidal","linear"),
    n_missing = 1,
    min_distance = NULL) {

  if (is.null(min_distance)) {
    min_distance <- round(1 / 5 * grid_points)
  }

  # Create grid points
  t <- seq(0, 10, length.out = grid_points) #seq(0, 1, length.out = grid_points)

  # Initialize storage
  curves <- vector("list", n)
  noisy_curves <- vector("list", n)
  # stochastic_components <- vector("list", n)
  if (!univariate) {
    curves_2 <- vector("list", n)
    noisy_curves_2 <- vector("list", n)
    # stochastic_components_2 <- vector("list", n)
  }
  nu <- numeric(n)

  # Helper function to generate curve with given parameters
  generate_curve_1 <- function(t, noise_sd, center = TRUE) {

      # Generate random coefficients for each subject
      u_i1 <- rnorm(1, mean = 0, sd = 5)  # u_i1 ~ N(0, 25)
      u_i2 <- rnorm(1, mean = 0, sd = 0.2) # u_i2 ~ N(0, 0.04)

      v_i1k <- rnorm(10, mean = 0, sd = 1) # v_i1k ~ N(0, 1)
      v_i2k <- rnorm(10, mean = 0, sd = 1) # v_i2k ~ N(0, 1)

      true_curve=rep(0,length(t))

      # For each grid point, compute X_i(t_g)
      for (g in 1:length(t)) {
        t_g <- t[g]

        # Base components
        base <- u_i1 + u_i2 * t_g

        # Sum components with sine and cosine terms
        sum_component <- 0
        for (k in 1:10) {
          sum_component <- sum_component +
            v_i1k[k] * sin((2*pi*k/10) * t_g) +
            v_i2k[k] * cos((2*pi*k/10) * t_g)
        }

        true_curve[g] <- base + sum_component
      }

      if (center == TRUE) {
        sol=list(
          curve_true = true_curve - mean(true_curve),
          curve_noisy = true_curve - mean(true_curve) + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
        )
      }else{
        sol=list(
          curve_true = true_curve,
          curve_noisy = true_curve + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
        )
      }

    sol
  }

  generate_curve_2 <- function(t, noise_sd, center = TRUE) {

    u_i1 <- rnorm(1, mean = 2, sd = 3)  # Different distribution than X1
    u_i2 <- rnorm(1, mean = -0.1, sd = 0.15) # Different slope distribution

    # Use different random seeds for v coefficients
    v_i1k <- rnorm(10, mean = -0.5, sd = 0.8) # Different than X1
    v_i2k <- rnorm(10, mean = 0.5, sd = 0.8)  # Different than X1

    true_curve=rep(0,length(t))

    # For each grid point, compute X2_i(t_g)
    for (g in 1:length(t)) {
      t_g <- t[g]

      # Different base components
      base <- u_i1 + u_i2 * t_g

      # Sum components with different frequencies and phase shifts
      sum_component <- 0
      for (k in 1:10) {
        # Use different frequencies and phase shifts for low correlation
        sum_component <- sum_component +
          v_i1k[k] * sin((2*pi*(k+0.5)/10) * t_g + pi/4) + # Phase shift by pi/4
          v_i2k[k] * cos((2*pi*(k+0.5)/10) * t_g - pi/4)   # Phase shift by -pi/4
      }

      true_curve[g] <- base + sum_component
    }

    if (center == TRUE) {
      sol=list(
        curve_true = true_curve - mean(true_curve),
        curve_noisy = true_curve - mean(true_curve) + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
      )
    }else{
      sol=list(
        curve_true = true_curve,
        curve_noisy = true_curve + stats::rnorm(length(t), 0, noise_sd * stats::sd(true_curve))
      )
    }
    sol
  }

  beta1 <- function(t) 0.05 * sin(pi*t/5)  # Scaled by 0.05
  beta2 <- function(t) 0.05 * (t/2.5)^2    # Scaled by 0.05


  # Generate beta coefficients with your original functions
  beta <- beta1(t)



  # Generate data and compute response
  for (i in 1:n) {

    # Generate curve
    curve_data <- generate_curve_1(
      t,
      noise_sd,
      center
    )

    curves[[i]] <- curve_data$curve_true
    noisy_curves[[i]] <- curve_data$curve_noisy

    if (univariate) {

      if (linear_predictor=="trapezoidal") {
        integrand <- curves[[i]] * beta
        nu[i] <- mu + (0.5 * sum(diff(t) * (integrand[-1] + integrand[-length(integrand)])))
      }else if(linear_predictor=="linear") {
        nu[i] <- mu + ((curves[[i]] %*% beta) / grid_points)
      }else if(linear_predictor=="rectangular"){
        integral1_approx = mean(curves[[i]] * beta) * (max(t) - min(t))
        nu[i] = mu + integral1_approx
      }

    }else{

      # Generate curve
      curve_data_2 <- generate_curve_2(
        t,
        noise_sd,
        center
      )

      curves_2[[i]] <- curve_data_2$curve_true
      noisy_curves_2[[i]] <- curve_data_2$curve_noisy

      beta_2 <- beta2(t)
      # beta_2_scaled <- scale_to_range(beta_2)


      integrand <- curves[[i]] * beta
      integrand_2 <- curves_2[[i]] * beta_2

      if (linear_predictor=="trapezoidal") {
        nu[i] <- mu + ( 0.5 * sum(diff(t) * (integrand[-1] + integrand[-length(integrand)])) +
                          0.5 * sum(diff(t) * (integrand_2[-1] + integrand_2[-length(integrand_2)])))
      }else if(linear_predictor=="linear"){
        nu[i] <- mu + (((curves[[i]] %*% beta) + (curves_2[[i]] %*% beta_2)) / grid_points)
      }else if(linear_predictor=="rectangular"){
        integral1_approx = mean(curves[[i]] * beta) * (max(t) - min(t))
        integral2_approx = mean(curves_2[[i]] * beta_2) * (max(t) - min(t))

       nu[i] = mu + integral1_approx + integral2_approx
      }

    }
  }

  if (response_type=="gaussian") {
    # Generate response with desired R-squared
    sd_e <- (1 / rsq - 1) * stats::sd(nu)
    # print(sd_e)
    response <- nu + stats::rnorm(n, 0, (sd_e))
  }else{
    response <- stats::rbinom(n, 1, (exp(nu) / (1 + exp(nu))))
    # if (sum(response)>=65 || sum(response)<=35) {
    #   zeros <- round(n/2)
    #   response <- c(rep(0,zeros),rep(1,n-zeros))
    # }
  }

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

  if (!univariate) {
    noisy_curves_miss_2 <- add_miss1d(
      noisy_curves_2,
      n_missing = n_missing,
      min_distance = min_distance
    )

    miss_points_2 <- noisy_curves_miss_2[["miss_points"]]
    missing_points_2 <- noisy_curves_miss_2[["missing_points"]]
    X_miss_2 <- matrix(
      unlist(noisy_curves_miss_2[["X_miss"]]),
      nrow = length(noisy_curves_miss_2[["X_miss"]]),
      ncol = length(noisy_curves_miss_2[["X_miss"]][[1]]),
      byrow = TRUE
    )
  }



  # Return results
  return(if (univariate) {
    list(
      curves = matrix(unlist(curves), nrow = n, byrow = TRUE),
      noisy_curves = matrix(unlist(noisy_curves), nrow = n, byrow = TRUE),
      noisy_curves_miss = X_miss,
      miss_points = miss_points,
      missing_points = missing_points,
      response = response,
      grid = t,
      beta = beta
    )

  }else{
    list(
      curves_1 = matrix(unlist(curves), nrow = n, byrow = TRUE),
      noisy_curves_1 = matrix(unlist(noisy_curves), nrow = n, byrow = TRUE),
      noisy_curves_miss_1 = X_miss,
      miss_points_1 = miss_points,
      missing_points_1 = missing_points,
      curves_2 = matrix(unlist(curves_2), nrow = n, byrow = TRUE),
      noisy_curves_2 = matrix(unlist(noisy_curves_2), nrow = n, byrow = TRUE),
      noisy_curves_miss_2 = X_miss_2,
      miss_points_2 = miss_points_2,
      missing_points_2 = missing_points_2,
      beta = beta,
      beta_2 = beta_2,
      grid = t,
      response = response
    )

  }
  )
}


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
    if (n_missing >= 1) {
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
#' @param n Number of samples to generate.
#' @param grid_x Number of points in x-axis grid. Default is 20.
#' @param grid_y Number of points in y-axis grid. Default is 20.
#' @param noise_sd Standard deviation of measurement noise. Default is 0.015.
#' @param rsq Desired R-squared value for the response. Default is 0.95.
#' @param beta_type Type of coefficient surface ("saddle" or "exp"). Default is "saddle".
#' @param response_type Type of the response variable ("gaussian" or "binomial"). Default is "gaussian".
#' @param a1 Optional fixed value for first stochastic component. If provided, a2 must also be provided.
#' @param a2 Optional fixed value for second stochastic component. If provided, a1 must also be provided.
#' @param sub_response Number of intervals for Simpson integration. Default is 50.
#' @param n_missing Number of holes in every curve.
#' @param min_distance_x Length of the holes in the x axis.
#' @param min_distance_y Length of the holes in the y axis.
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
#' @examples
#' # Generate basic 2D functional data with default parameters
#' data <- data_generator_po_2d(n = 2)
#'
#' # Generate data with custom grid size and Gaussian response
#' data <- data_generator_po_2d(n = 2, grid_x = 30, grid_y = 30, response_type = "gaussian")
#'
#' # Generate data with binomial response and saddle-shaped coefficient surface
#' data <- data_generator_po_2d(n = 2, response_type = "binomial", beta_type = "saddle")
#'
#' # Generate data with fixed stochastic components
#' data <- data_generator_po_2d(n = 2, a1 = 0.1, a2 = -0.2)
#'
#' # Introduce missing data with holes along curves
#' data <- data_generator_po_2d(n = 2, n_missing = 3, min_distance_x = 5, min_distance_y = 5)
#'
#' @export
data_generator_po_2d <- function(
    n = 20,
    grid_x = 20,
    grid_y = 20,
    noise_sd = 0.015,
    rsq = 0.95,
    intercept = 0.1,
    beta_type = c("saddle", "exp", "smooth", "sinusoidal", "peaks"),
    response_type = c("gaussian", "binomial"),
    linear_predictor = c("integral","linear"),
    a1 = NULL,
    a2 = NULL,
    sub_response = 50,
    n_missing = 1,
    min_distance_x = NULL,
    min_distance_y = NULL) {
  beta_type <- match.arg(beta_type)
  response_type <- match.arg(response_type)

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
    min_distance_x <- round(1 / 5 * grid_x)
  }

  if (is.null(min_distance_y)) {
    min_distance_y <- round(1 / 5 * grid_y)
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


  # Prepare integration grid
  n_x <- m_y <- 2 * sub_response

  x_fine <- seq(min(x), max(x), length.out = n_x + 1)
  y_fine <- seq(min(y), max(y), length.out = m_y + 1)

  h_x <- (max(x_fine) - min(x_fine)) / n_x
  h_y <- (max(y_fine) - min(y_fine)) / m_y


  # Generate beta surface on fine grid
  if(linear_predictor == "integral"){
    beta_fine <- if (beta_type == "saddle") {
      generate_saddle_surface(x_fine, y_fine)
    }else if(beta_type == "exp") {
      generate_exp_surface(x_fine, y_fine)
    }else if(beta_type == "smooth"){
      generate_smooth_surface(x_fine, y_fine)
    }else if(beta_type == "peaks"){
     generate_multipeak_surface(x_fine, y_fine)
    }
    else{
      generate_sinusoidal_surface(x_fine, y_fine)
    }}else{
      beta_fine <- if (beta_type == "saddle") {
        generate_saddle_surface(x, y)
      }else if(beta_type == "exp") {
        generate_exp_surface(x, y)
      }else if(beta_type == "smooth"){
        generate_smooth_surface(x, y)
      }else if(beta_type == "peaks"){
        generate_multipeak_surface(x, y)
      }else{
        generate_sinusoidal_surface(x, y)
      }}

  # Generate beta on original grid for output
  beta_surface <- if (beta_type == "saddle") {
    generate_saddle_surface(x, y)
  }else if(beta_type == "exp") {
    generate_exp_surface(x, y)
  }else if(beta_type == "smooth"){
    generate_smooth_surface(x, y)
  }else if(beta_type == "peaks"){
    generate_multipeak_surface(x, y)
  }else{
    generate_sinusoidal_surface(x, y)
  }


  # Generate data and compute response
  for (i in 1:n) {
    # Use fixed values if provided, otherwise generate random components
    if (!is.null(a1)) {
      stochastic_components[i, ] <- c(a1, a2)
    } else {
      stochastic_components[i, ] <- c(
        stats::rnorm(1, 0, 1), # a1
        stats::rnorm(1, 0, 1) # a2
      )
    }

    # Generate surface
    surface_data <- generate_surface(
      x, y,
      stochastic_components[i, 1],
      stochastic_components[i, 2],
      noise_sd
    )

    surfaces[[i]] <- surface_data$DATA_T
    noisy_surfaces[[i]] <- surface_data$DATA_N




    # Generate finer surface for integration
    if(linear_predictor == "integral"){

      surface_fine <- generate_surface(
      x_fine, y_fine,
      stochastic_components[i, 1],
      stochastic_components[i, 2],
      0
    )$DATA_T
    }

    if(linear_predictor == "integral"){

      integrand = surfaces[[i]] %*% beta_surface

      nu[i] <- double_integral(integrand, x, y)

      # Get integration weights
      # W_delta <- setup_simpson_weights(n_x, m_y, h_x, h_y)

      # Compute double integral using Simpson's rule
      # nu[i] <- as.double(t(as.vector(surface_fine)) %*%
      #                      diag(W_delta) %*%
      #                      as.vector(beta_fine))
    }else{

      nu[i] <- mean(surfaces[[i]] * beta_surface)

      # nu[i] <- as.double(t(as.vector(surfaces[[i]])) %*%
      #                      as.vector(beta_surface))
    }

  }

  if (response_type == "gaussian") {
    # Generate response with desired R-squared
    var_e <- (1 / rsq - 1) * stats::var(nu)
    response <- intercept + nu + stats::rnorm(n, 0, sqrt(var_e))
  } else {

    # stats::rbinom(n, 1, (exp(nu) / (1 + exp(nu))))

    aux=adjust_proportion(intercept, nu, max_iter = 150, tolerance = 0.001)

    response = aux$y
    intercept = aux$intercept

  }

  noisy_surfaces_miss <- add_miss2(
    noisy_surfaces,
    n_missing,
    min_distance_x,
    min_distance_y
  )

  list(
    surfaces = surfaces,
    noisy_surfaces = noisy_surfaces,
    noisy_surfaces_miss = noisy_surfaces_miss[[1]],
    miss_points = noisy_surfaces_miss[[2]],
    missing_points = noisy_surfaces_miss[[3]],
    response = response,
    nu = nu,
    intercept = intercept,
    points_x = x,
    points_y = y,
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

#' Generate smooth coefficient surface
#'
#' @noRd
generate_smooth_surface <- function(x, y) {
  surface <- matrix(nrow = length(x), ncol = length(y))

  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      surface[i, j] <- 0.5* (x[i])^2 +  0.5 * (y[j])^2
    }
  }
  surface
}

#' Generate sinusoidal coefficient surface with Gaussian decay
#'
#' @noRd
generate_sinusoidal_surface <- function(x, y) {
  surface <- matrix(nrow = length(x), ncol = length(y))
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      surface[i, j] <- 0.8 * sin(2 * pi * x[i]) * cos(2 * pi * y[j]) *
        exp(-2 * (x[i] - 0.5)^2) * exp(-2 * (y[j] - 0.5)^2)
    }
  }
  surface
}

#' Generate multi-peak coefficient surface
#'
#' @noRd
generate_multipeak_surface <- function(x, y) {
  surface <- matrix(nrow = length(x), ncol = length(y))
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      surface[i, j] <- 0.5 * sin(4 * pi * x[i]) * cos(4 * pi * y[j]) +
        0.3 * exp(-5 * ((x[i] - 0.3)^2 + (y[j] - 0.7)^2))
    }
  }
  surface
}

#' Generate simple test coefficient surface
#'
#' @noRd
generate_test_surface <- function(x, y) {
  surface <- matrix(nrow = length(x), ncol = length(y))
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      # Simple quadratic surface
      surface[i, j] <- 5 * ((x[i] - 0.5)^2 + (y[j] - 0.5)^2) - 2
      # surface[i, j] <- 2 * (x[i] - 0.5)^2 + 2 * (y[j] - 0.5)^2 - 0.5
    }
  }
  surface
}
#' Helper function to generate surface X(t_1,t_2) with given parameters
#'
#' @noRd
generate_surface <- function(x, y, a1, a2, noise_sd) {
  true_surface <- matrix(nrow = length(x), ncol = length(y))
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      true_surface[i, j] <- a1 * cos(2 * pi * x[i]) +
        a2 * sin(2 * pi * y[j]) + 1 #+ 2 * y[j] - 0.5 * x[i]*y[j]
    }
  }

  list(
    DATA_T = true_surface,
    DATA_N = true_surface + matrix(
      stats::rnorm(length(x) * length(y), 0, noise_sd*sd(as.vector(true_surface))),
      length(x),
      length(y)
    )
  )
}

#' Generate 2D functional predictors using B-spline basis expansion
#'
#' @noRd
generate_bspline_surfaces <- function(n_obs, x, y, n_basis_x = 8, n_basis_y = 8, coef_sd = 0.3, noise_sd) {

  # Create tensor product B-spline basis for 2D surfaces
  basis_x <- fda::create.bspline.basis(c(0, 1), n_basis_x)
  basis_y <- fda::create.bspline.basis(c(0, 1), n_basis_y)

  # Generate coefficients for each observation's surface
  # Each observation has a 2D surface X_i(t_1,t_2)
  X_coefs <- array(stats::rnorm(n_obs * n_basis_x * n_basis_y, 0, coef_sd),
                   dim = c(n_obs, n_basis_x, n_basis_y))

  # Evaluate surfaces on the grid manually using tensor products
  X_surfaces <- array(0, dim = c(n_obs, length(x), length(y)))
  X_noisy_surfaces <- array(0, dim = c(n_obs, length(x), length(y)))

  # Evaluate basis functions on grids
  basis_x_vals <- fda::eval.basis(x, basis_x)  # n_s x n_basis_x
  basis_y_vals <- fda::eval.basis(y, basis_y)  # n_t x n_basis_y

  for(i in 1:n_obs) {
    # Compute tensor product: X_i(t_1,t_2) = sum_j sum_k c_ijk * phi_j(t_1) * psi_k(t_2)
    for(j in 1:n_basis_x) {
      for(k in 1:n_basis_y) {
        # Add contribution of basis function (j,k) with coefficient c_ijk
        X_surfaces[i,,] <- X_surfaces[i,,] +
          X_coefs[i,j,k] * outer(basis_x_vals[,j], basis_y_vals[,k])
      }
    }

    # Add noise to create observed surface
    X_noisy_surfaces[i,,] <- X_surfaces[i,,] +
      matrix(stats::rnorm(length(x) * length(y), 0,
                          noise_sd * stats::sd(as.vector(X_surfaces[i,,]))),
             length(x), length(y))
  }

  list(
    surfaces = X_surfaces,
    noisy_surfaces = X_noisy_surfaces,
    coefficients = X_coefs,
    basis_x = basis_x,
    basis_y = basis_y
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
  W_x <- (h_x / 3) * simp_w_x

  # Create full weight matrix
  for (j in 1:(n_y + 1)) {
    start_idx <- ((n_x + 1) * (j - 1) + 1)
    end_idx <- ((n_x + 1) * j)

    if (j == 1 || j == (n_y + 1)) {
      W_delta[start_idx:end_idx] <- (h_y / 3) * W_x
    } else if (j %% 2 == 0) {
      W_delta[start_idx:end_idx] <- (4 * h_y / 3) * W_x
    } else {
      W_delta[start_idx:end_idx] <- (2 * h_y / 3) * W_x
    }
  }

  W_delta
}

#' Setup 2D Trapeizodal integration
#'
#' @noRd
double_integral <- function(integrand, x, y) {
  # Get dimensions
  nx <- nrow(integrand)
  ny <- ncol(integrand)

  dx <- diff(x)[1]  # spacing in x direction
  dy <- diff(y)[1]  # spacing in y direction

  # Corner points (weight = 1/4)
  corner_sum <- sum(integrand[1,1], integrand[1,ny],
                   integrand[nx,1], integrand[nx,ny],na.rm=TRUE) / 4

  # Edge points (weight = 1/2)
  edge_sum <- (sum(integrand[1, 2:(ny-1)],na.rm=TRUE) +    # top edge
                 sum(integrand[nx, 2:(ny-1)],na.rm=TRUE) +    # bottom edge
                 sum(integrand[2:(nx-1), 1],na.rm=TRUE) +     # left edge
                 sum(integrand[2:(nx-1), ny],na.rm=TRUE)) / 2 # right edge

  # Interior points (weight = 1)
  interior_sum <- sum(integrand[2:(nx-1), 2:(ny-1)],na.rm=TRUE)

  # Combine all parts and multiply by grid spacing
  result <- dx * dy * (corner_sum + edge_sum + interior_sum)

  return(result)
}

#' Iteratively adjust intercept to achieve target proportion in binomial simulation
#'
#' This function uses an iterative approach to find the appropriate intercept value
#' that produces a desired proportion of 1s in binomial response simulation. It works
#' by adjusting the intercept on the log-odds scale using adaptive damping to prevent
#' overshooting due to the nonlinear logistic transformation.
#'
#' @param target_prop Desired proportion of 1s in the response. Must be between 0 and 1 (exclusive).
#' @param functional_effects Numeric vector of functional effects (e.g., from 2D integration
#'   of surfaces). These represent the variability around the baseline intercept.
#' @param max_iter Maximum number of iterations for adjustment. Default is 15.
#' @param tolerance Convergence tolerance for the difference between target and achieved
#'   proportion. Default is 0.03.
#' @param verbose Logical indicating whether to print iteration progress. Default is FALSE.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Starting with an initial intercept based on \code{qlogis(target_prop)}
#'   \item Computing probabilities using \code{plogis(intercept + functional_effects)}
#'   \item Generating binary outcomes using \code{rbinom()}
#'   \item Adjusting the intercept based on the error between target and achieved proportions
#'   \item Using adaptive damping (0.3 for large errors, 0.5 for medium, 0.7 for small)
#' }
#'
#' The adjustment formula is:
#' \code{adjustment = (qlogis(target_prop) - qlogis(achieved_prop)) * damping_factor}
#'
#' This approach works in log-odds space to prevent probabilities from exceeding [0,1]
#' bounds and provides robust control over response proportions in functional regression
#' simulation studies.
#'
#' @return A list containing:
#' \itemize{
#'   \item y: Binary response vector of length equal to \code{functional_effects}
#'   \item intercept: Final adjusted intercept value
#'   \item final_prop: Achieved proportion of 1s in the response
#'   \item iterations: Number of iterations used
#'   \item converged: Logical indicating whether convergence was achieved
#'   \item final_error: Final absolute error between target and achieved proportion
#' }
#'
#' @examples
#' # Basic usage with simulated functional effects
#' set.seed(123)
#' effects <- rnorm(100, mean = 0, sd = 0.5)
#' result <- adjust_proportion(target_prop = 0.3, functional_effects = effects)
#' cat("Achieved proportion:", result$final_prop, "\n")
#' cat("Converged in", result$iterations, "iterations\n")
#'
#' # Usage with 2D functional regression effects
#' # Assuming you have computed 2D integral effects from surfaces
#' # integral_effects <- compute_2d_integrals(surfaces, beta_surface)
#' # result <- adjust_proportion(0.4, integral_effects, tolerance = 0.01)
#'
#' # Check for convergence issues
#' if (!result$converged) {
#'   warning("Adjustment did not converge. Final error: ", result$final_error)
#' }
#'
#' @seealso \code{\link{data_generator_po_2d}} for using this function in 2D functional
#'   data simulation.
adjust_proportion <- function(target_prop, functional_effects,
                              max_iter = 15, tolerance = 0.03, verbose = FALSE) {

  # Input validation
  if(any(is.na(functional_effects)) || any(is.infinite(functional_effects))) {
    stop("functional_effects contains NA or infinite values", call. = FALSE)
  }

  if(target_prop <= 0 || target_prop >= 1) {
    stop("target_prop must be between 0 and 1 (exclusive)", call. = FALSE)
  }

  if(max_iter < 1) {
    stop("max_iter must be at least 1", call. = FALSE)
  }

  if(tolerance <= 0) {
    stop("tolerance must be positive", call. = FALSE)
  }

  intercept <- qlogis(target_prop)

  # Check if initial intercept is reasonable
  if(is.infinite(intercept)) {
    stop("target_prop too close to 0 or 1, causing infinite initial intercept", call. = FALSE)
  }

  for(iter in 1:max_iter) {
    linear_pred <- intercept + functional_effects

    # Check for extreme linear predictors
    if(any(is.na(linear_pred)) || any(is.infinite(linear_pred))) {
      if(verbose) cat("Warning: NA or infinite linear predictors at iteration", iter, "\n")
      # Try to recover by reducing intercept magnitude
      intercept <- intercept * 0.5
      next
    }

    # Check for extreme linear predictors that would cause numerical issues
    if(max(abs(linear_pred)) > 20) {
      if(verbose) cat("Warning: Very large linear predictors (max =", max(abs(linear_pred)), ") at iteration", iter, "\n")
      # Clip extreme values
      linear_pred <- pmax(pmin(linear_pred, 20), -20)
    }

    probabilities <- plogis(linear_pred)

    # Check probabilities for issues
    if(any(is.na(probabilities))) {
      if(verbose) cat("Warning: NA probabilities at iteration", iter, "\n")
      intercept <- intercept * 0.8
      next
    }

    # Generate binary outcomes
    y <- rbinom(length(functional_effects), 1, probabilities)

    # Check for NA in y
    if(any(is.na(y))) {
      if(verbose) cat("Warning: NA in binary outcomes at iteration", iter, "\n")
      intercept <- intercept * 0.8
      next
    }

    achieved_prop <- mean(y)

    # Check if achieved_prop is valid
    if(is.na(achieved_prop)) {
      if(verbose) cat("Warning: NA achieved proportion at iteration", iter, "\n")
      intercept <- intercept * 0.8
      next
    }

    # Handle edge cases BEFORE computing error
    # This prevents qlogis() from producing infinite values
    if(achieved_prop == 0) {
      achieved_prop <- 0.5 / length(functional_effects)  # Small positive value
    } else if(achieved_prop == 1) {
      achieved_prop <- 1 - 0.5 / length(functional_effects)  # Small value less than 1
    }

    error <- abs(achieved_prop - target_prop)

    # Check if error is valid
    if(is.na(error)) {
      if(verbose) cat("Warning: NA error at iteration", iter, "\n")
      break
    }

    # Adaptive damping based on error size
    if(error > 0.1) {
      damp <- 0.3  # Very conservative for large errors
    } else if(error > 0.05) {
      damp <- 0.5  # Moderate for medium errors
    } else {
      damp <- 0.7  # More aggressive for small errors
    }

    if(verbose) cat("Iteration", iter, ": achieved_prop =", round(achieved_prop, 4),
                    ", error =", round(error, 4), ", damping =", damp, "\n")

    # Check convergence
    if(error <= tolerance) {
      if(verbose) cat("Converged after", iter, "iterations\n")
      break
    }

    # Normal adjustment (no more edge case handling needed)
    adjustment <- (qlogis(target_prop) - qlogis(achieved_prop)) * damp
    intercept <- intercept + adjustment

    # Check if new intercept is reasonable
    if(is.infinite(adjustment) || is.na(adjustment)) {
      if(verbose) cat("Warning: Invalid adjustment at iteration", iter, "\n")
      # Use a small fixed adjustment instead
      if(achieved_prop < target_prop) {
        intercept <- intercept + 0.1
      } else {
        intercept <- intercept - 0.1
      }
    }

    # Prevent intercept from becoming too extreme
    if(abs(intercept) > 10) {
      if(verbose) cat("Warning: Intercept becoming too extreme (", intercept, "), clipping\n")
      intercept <- sign(intercept) * 10
    }
  }

  # Final check
  if(iter >= max_iter) {
    warning("Maximum iterations reached without convergence. Final error: ", round(error, 4))
  }

  return(list(
    y = y,
    intercept = intercept,
    final_prop = achieved_prop,
    iterations = iter,
    converged = error <= tolerance,
    final_error = error
  ))
}


#' Generate high-signal 2D functional data with proper interface
#'
#' Creates synthetic 2D functional data using the working high-signal approach
#' but with the same interface as data_generator_po_2d for compatibility
#'
#' @param n Number of samples to generate
#' @param grid_x Number of points in x-axis grid
#' @param grid_y Number of points in y-axis grid
#' @param intercept Target proportion for binomial or intercept for Gaussian
#' @param noise_sd Standard deviation of measurement noise
#' @param beta_type Type of coefficient surface (kept for compatibility)
#' @param response_type Type of response ("binomial" or "gaussian")
#' @param signal_strength Strength of discriminative signal (default 2.5)
#' @param n_missing Not used (kept for compatibility)
#' @param min_distance_x Not used (kept for compatibility)
#' @param min_distance_y Not used (kept for compatibility)
#' @param linear_predictor Not used (kept for compatibility)
#' @param sub_response Not used (kept for compatibility)
#'
#' @return A list containing simulated data with same structure as data_generator_po_2d

#' @export
data_generator_high_signal <- function(n = 100,
                                       grid_x = 20,
                                       grid_y = 20,
                                       intercept = 0.6,
                                       noise_sd = 0.25,
                                       response_type = "binomial",
                                       signal_strength = 2.5,
                                       n_missing = 0,
                                       min_distance_x = NULL,
                                       min_distance_y = NULL){

  # Create coordinate grids
  x <- seq(0, 1, length.out = grid_x)
  y <- seq(0, 1, length.out = grid_y)

  if(response_type == "binomial") {
    target_prop <- intercept
  } else {
    target_prop <- 0.5  # Not used for Gaussian
  }

  cat("=== HIGH-SIGNAL DATA GENERATION ===\n")
  cat("Response type:", response_type, "\n")
  cat("Sample size:", n, "\n")
  cat("Grid size:", grid_x, "x", grid_y, "\n")
  if(response_type == "binomial") {
    cat("Target proportion:", target_prop, "\n")
  }

  # Step 1: Create discriminative coefficient surface
  beta_surface <- matrix(nrow = length(x), ncol = length(y))
  for(i in seq_along(x)) {
    for(j in seq_along(y)) {
      # Strong discriminative pattern
      beta_surface[i, j] <- signal_strength * (
        sin(2*pi*x[i]) * cos(2*pi*y[j]) +
          0.8 * (x[i] - 0.5) * (y[j] - 0.5)
      )
    }
  }

  cat("Beta surface range:", range(beta_surface), "\n")

  # Step 2: Generate functional predictor surfaces
  surfaces <- vector("list", n)
  noisy_surfaces <- vector("list", n)
  stochastic_components <- matrix(nrow = n, ncol = 2,
                                  dimnames = list(NULL, c("a1", "a2")))
  functional_effects <- numeric(n)

  group1_size <- round(n * target_prop)
  group_assignment <- c(rep(1, group1_size), rep(2, n - group1_size))
  group_assignment <- sample(group_assignment)

  for(i in 1:n) {
    # Generate coefficients based on group membership
    if(group_assignment[i] == 1) {
      # Group 1: positive correlation with beta surface
      a1 <- rnorm(1, 1.5, 0.3)
      a2 <- rnorm(1, 1.2, 0.3)
    } else {
      # Group 2: negative correlation with beta surface
      a1 <- rnorm(1, -1.2, 0.3)
      a2 <- rnorm(1, -1.0, 0.3)
    }

    stochastic_components[i, ] <- c(a1, a2)

    # Generate surface using matching basis functions
    true_surface <- matrix(nrow = length(x), ncol = length(y))
    for(ii in seq_along(x)) {
      for(jj in seq_along(y)) {
        true_surface[ii, jj] <-
          a1 * sin(2*pi*x[ii]) +                    # Matches beta's sin component
          a2 * cos(2*pi*y[jj]) +                    # Matches beta's cos component
          0.5 * (x[ii] - 0.5) * (y[jj] - 0.5) +    # Matches beta's interaction
          2.0                                        # Baseline level

        ## THIS ONE IS TO USE WITH THE double_integral() FUNCTION
        # true_surface[ii, jj] <-
        #   a1 * (sin(2*pi*x[ii]) + 1) +                    # Shift to positive
        #   a2 * (cos(2*pi*y[jj]) + 1) +                    # Shift to positive
        #   0.5 * (x[ii] - 0.5) * (y[jj] - 0.5) +
        #   2.0
      }
    }

    surfaces[[i]] <- true_surface

    # Add noise proportional to surface variation
    surface_sd <- sd(as.vector(true_surface))
    noise_matrix <- matrix(rnorm(length(x) * length(y), 0, noise_sd * surface_sd),
                           length(x), length(y))
    noisy_surfaces[[i]] <- true_surface + noise_matrix

    # Compute functional effect (integral)
    integrand <- true_surface * beta_surface
    functional_effects[i] <- mean(integrand) * (max(x) - min(x)) * (max(y) - min(y))
    # functional_effects[i] <- double_integral(integrand, x, y)
  }

  cat("Functional effects range:", range(functional_effects), "\n")
  cat("Functional effects std:", sd(functional_effects), "\n")

  # Step 3: Generate response based on type
  if(response_type == "binomial") {
    # Binomial response using functional effects
    # scaled_effects <- scale(functional_effects)[,1] * 1.5
    scaled_effects <- functional_effects

    probabilities <- plogis(qlogis(target_prop) + scaled_effects)
    response <- rbinom(n, 1, probabilities)

    cat("Response proportion:", mean(response), "\n")
    final_intercept <- qlogis(target_prop)

  } else {
    # Gaussian response
    rsq <- 0.95
    var_e <- (1 / rsq - 1) * var(functional_effects)
    response <- intercept + functional_effects + rnorm(n, 0, sqrt(var_e))

    cat("Response range:", range(response), "\n")
    cat("Response mean:", mean(response), "\n")
    final_intercept <- intercept
  }

  noisy_surfaces_miss <- add_miss2(
    noisy_surfaces,
    n_missing,
    min_distance_x,
    min_distance_y
  )


  # Step 4: Create data structure matching data_generator_po_2d output
  simulation_data <- list(
    # Core surfaces (same format as original)
    surfaces = surfaces,
    noisy_surfaces = noisy_surfaces,
    noisy_surfaces_miss = noisy_surfaces_miss[[1]],


    # Missing data structures (empty for now)
    miss_points = noisy_surfaces_miss[[2]],
    missing_points = noisy_surfaces_miss[[3]],

    # Response and parameters
    response = response,
    intercept = final_intercept,

    # Grid information (matching original naming)
    points_x = x,
    points_y = y,

    # True parameters
    beta = beta_surface,
    stochastic_components = stochastic_components,

    # Additional information
    functional_effects = functional_effects,
    group_assignment = group_assignment,
    response_type = response_type,
    signal_strength = signal_strength
  )

  # Step 5: Diagnostics
  if(response_type == "binomial") {
    correlation <- cor(functional_effects, response)
    simple_model <- glm(response ~ functional_effects, family = binomial())
    simple_pred <- predict(simple_model, type = "response")

    simple_auc <- tryCatch({
      if(requireNamespace("pROC", quietly = TRUE)) {
        as.numeric(pROC::auc(response, simple_pred, quiet = TRUE))
      } else {
        NA
      }
    }, error = function(e) NA)

    cat("Correlation (functional effects vs response):", round(correlation, 3), "\n")
    cat("Simple model AUC:", round(simple_auc, 3), "\n")

    if(correlation > 0.3 && !is.na(simple_auc) && simple_auc > 0.7) {
      cat("SUCCESS: Strong discriminative signal achieved!\n")
    } else {
      cat("WARNING: Signal may be insufficient\n")
    }

    simulation_data$diagnostics <- list(
      correlation = correlation,
      simple_auc = simple_auc
    )
  } else {
    # Gaussian diagnostics
    correlation <- cor(functional_effects, response)
    cat("Correlation (functional effects vs response):", round(correlation, 3), "\n")

    simulation_data$diagnostics <- list(
      correlation = correlation,
      r_squared = correlation^2
    )
  }

  return(simulation_data)
}
