#' Data generator function for the variable domain case.
#'
#' This function is for internal use.
#'
#' @param N Number of subjects.
#' @param J Number of maximum observations per subject.
#' @param nsims Number of simulations per the simulation study.
#' @param aligned If the data that will be generated is aligned or not.
#' @param multivariate If TRUE, the data is generated with 2 variables.
#' @param beta_index Index for the beta.
#' @param Rsq Variance of the model.
#' @param use_x If the data is generated with x.
#' @param use_f If the data is generated with f.
#' @param seed Seed for reproducibility.
#'
#' @return Example data.
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
    use_f = FALSE,
    seed = NULL
) {
  if (!(beta_index %in% c(1, 2))) {
    stop("'beta_index' could only be 1 or 2",call. = FALSE)
  }

  for (iter in 1:nsims) {
    if (!is.null(seed)) {
      set.seed(seed + iter)
    } else {
      set.seed(iter)
    }

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

    # if (max(M) > J) {
    #   M[which(M > J)] <- J
    # }
    #
    # if (min(M) <= 10) {
    #   M[which(M <= 10)] <- 10
    # }

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

      aux <- var(B, na.rm = TRUE)

      X_s[i, ] <- B
      X_se[i, ] <- B + stats::rnorm(maxM, 0, sqrt(aux / 8)) # WE ADD NOISE
      Y_s[i, ] <- B2
      Y_se[i, ] <- B2 + stats::rnorm(maxM, 0, sqrt(aux / 8)) # WE ADD NOISE
    }

    Beta <- array(dim = c(N, maxM, 4))
    nu <- rep(0, N)
    y <- rep(0, N)
    x1 <- stats::runif(N)
    x2 <- stats::rnorm(N)
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
        Beta[i, (M[i, 1]:M[i, 2]), 2] <- ((1 - (2 * M_diff[i] / T)) * (5 - 40 * ((t[(M[i, 1]:M[i, 2])] / M_diff[i]) - 0.5)^2)) / 10

        if (multivariate) {
          nu[i] <- sum(X_s[i, ] * Beta[i, ,beta_index], na.rm = TRUE) / (M_diff[i]) + sum(Y_s[i, ] * Beta[i, , 2], na.rm = TRUE) / (M_diff[i])
        } else {
          nu[i] <- sum(X_s[i, ] * Beta[i, ,beta_index], na.rm = TRUE) / (M_diff[i]) # NOT NOISY
        }

      }

    }

    nu <- if (use_f) nu + f1(x1) else nu
    var_e <- (1 / Rsq - 1) * stats::var(nu)
    y <- nu + stats::rnorm(N, sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL
    y <- if (use_x) y + x2 else y
  }


  data <- list(y = y)

  data[["X_s"]] <- X_s
  data[["X_se"]] <- X_se
  data[["Y_s"]] <- Y_s
  data[["Y_se"]] <- Y_se
  data[["x1"]] <- x1
  data[["x2"]] <- x2
  data[["Beta"]] <- Beta

  data
}
