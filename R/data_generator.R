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
#'
#' @return Example data.
#'
#' @noRd
data_generator_vd <- function(N = 100, J = 100, nsims = 1, Rsq = 0.95, aligned = TRUE, multivariate = FALSE, beta_index = 1, use_x = FALSE, use_f = FALSE) {
  if (!(beta_index %in% c(1, 2))) {
    stop("'beta_index' could only be 1 or 2",call. = FALSE)
  }

  for (iter in 1:nsims) {
    set.seed(42 + iter)

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
            v_i1 * sin(2 * pi * k * (1:M[i]) / 100) + v_i2 * cos(2 * pi * k * (1:M[i]) / 100)
        } else {
          temp[k, (M[i, 1]:M[i, 2])] <-
            v_i1 * sin(2 * pi * k * (M[i, 1]:M[i, 2]) / 100) + v_i2 * cos(2 * pi * k * (M[i, 1]:M[i, 2]) / 100)
        }
      }

      B <- apply(temp, 2, sum)


      B <- B + u1
      B2 <- B + stats::rnorm(1, sd = 0.02) + (t / 10)

      X_s[i, ] <- B
      X_se[i, ] <- B + stats::rnorm(maxM, 0, 1) # WE ADD NOISE
      Y_s[i, ] <- B2
      Y_se[i, ] <- B2 + stats::rnorm(maxM, 0, 1) # WE ADD NOISE
    }

    Beta <- array(dim = c(N, maxM, 4))
    nu <- rep(0, N)
    y <- rep(0, N)
    x1 <- stats::runif(N)
    f0 <- function(x) 2 * sin(pi * x)

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

    nu <- if (use_f) nu + f0(x1) else nu
    var_e <- (1 / Rsq - 1) * stats::var(nu)
    y <- nu + stats::rnorm(N, sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL
    y <- if (use_x) y + x1 else y
  }


  data <- data.frame(y = y)

  data[["X_s"]] <- X_s
  data[["X_se"]] <- X_se
  data[["Y_s"]] <- Y_s
  data[["Y_se"]] <- Y_se
  data[["x1"]] <- x1
  data[["Beta"]] <- Beta

  data
}


#' Data generator for the 1-dimensional partially observed functional data case
#'
#' @return Example data.
#'
#' @noRd
data_generator_po <- function() {
  data.gen.splines.uneven <- function(n, p, grid, knots, mean, err, coef) {
    ########################################################################################
    ##                                                                                    ##
    #       genera un campione di curve da una base spline con nodi irregolari
    #
    #       genera una muestra de curvas a partir de una base spline con nodos irregulares
    ##                                                                                    ##
    ########################################################################################

    #     n: number of subjects
    #     p: number of observation points
    #  grid: evaluation grid
    # knots: vector of knots placement
    #  mean: mean vector for the coefficients
    #   err: standard error of measurement error
    #  coef: variance vector for the coefficients (the covariance matrix is supposed to be diagonal)

    x <- matrix(nrow = n, ncol = p)
    x.sm <- matrix(nrow = n, ncol = p)

    basis <- fda::create.bspline.basis(rangeval = grid[c(1, p)], breaks = knots)

    for (i in 1:n) {
      coefs <- mvtnorm::rmvnorm(1, mean, coef * diag(rep(1, length(knots) + 2)))
      x.sm[i, ] <- as.vector(coefs %*% t(fda::eval.basis(grid, basis))) # smooth data (B-splines evaluated)
      x[i, ] <- as.vector(coefs %*% t(fda::eval.basis(grid, basis))) + stats::rnorm(p, sd = err) # noisy data (B-splines evaluated + an error)
    }

    out <- list(x, x.sm)
    names(out) <- c("data", "smooth.data")
    return(out)
  }
  data.agg.miss <- function(x) {
    ########################################################################################
    ##                                                                                    ##
    #       genera starting points ed ending points da una distribuzione uniforme          #
    ##
    #       genera una matriz con datos faltantes (las curvas partially observed)
    ##
    #       todas las curvas tienen el intervalo común timepoints
    ########################################################################################

    # x: data matrix, for example "rbind()" of the previous data generation

    p <- ncol(x)
    n <- nrow(x)
    init <- round(stats::runif(n - 1, 1, round(p / 3))) # generiamo i punti di inizio e fine
    fin <- round(stats::runif(n - 1, round(2 * p / 3 + 1), p))
    ep <- rbind(c(1, p), cbind(init, fin))

    for (i in 2:n) { # generiamo gli NAs
      aa <- ep[i, 1]
      bb <- ep[i, 2]
      x[i, 1:aa] <- NA
      x[i, bb:p] <- NA
    }

    obs <- matrix(nrow = n, ncol = 2)
    for (i in 1:n) {
      obs[i, ] <- c(min(which(!is.na(x[i, ]))), max(which(!is.na(x[i, ]))))
    }
    min <- max(obs[, 1])
    max <- min(obs[, 2])
    t <- min:max
    m <- length(t)

    domain <- rbind(cbind(1, p), cbind(ep[2:n, 1] + 1, ep[2:n, 2] - 1))

    out <- list(x, t, m, domain)
    names(out) <- c("x.miss", "timepoints", "interval.length", "domain")
    return(out)
  }
  data.agg.miss.betadecay <- function(x, rate) {
    ########################################################################################
    ##                                                                                    ##
    #         genera starting points ed ending points da una distribuzione beta            #
    ##                                                                                    ##
    ########################################################################################

    # x: data matrix, for example "rbind()" of the previous data generation

    p <- ncol(x)
    n <- nrow(x)
    init <- round((1 - stats::rbeta(n - 1, 1, rate)) * (p / 3)) # generiamo i punti di inizio e fine

    if (!is.null(which(init == 0))) {
      init[which(init == 0)] <- 1
    }

    fin <- round(2 * p / 3 + stats::rbeta(n - 1, 1, rate) * (p / 3))
    ep <- rbind(c(1, p), cbind(init, fin))

    for (i in 2:n) { # generiamo gli NAs
      aa <- ep[i, 1]
      bb <- ep[i, 2]
      x[i, 1:aa] <- NA
      x[i, bb:p] <- NA
    }

    obs <- matrix(nrow = n, ncol = 2)
    for (i in 1:n) {
      obs[i, ] <- c(min(which(!is.na(x[i, ]))), max(which(!is.na(x[i, ]))))
    }
    min <- max(obs[, 1])
    max <- min(obs[, 2])
    t <- min:max
    m <- length(t)

    domain <- rbind(cbind(1, p), cbind(ep[2:n, 1] + 1, ep[2:n, 2] - 1))

    out <- list(x, t, m, domain)
    names(out) <- c("x.miss", "timepoints", "interval.length", "domain")
    return(out)
  }

  p <- 150
  n1 <- 50
  n2 <- 50
  n <- n1 + n2
  grid <- seq(0, 1, length = p)
  err1 <- 0.1
  err2 <- 0.1
  unif.miss <- TRUE


  nknots <- 18
  par.knots <- nknots / 3

  knots <- c(
    sort(stats::runif(par.knots, grid[1], grid[p] / 3)),
    sort(stats::runif(par.knots, grid[p] / 3, 2 * grid[p] / 3)),
    sort(stats::runif(par.knots, 2 * grid[p] / 3, grid[p]))
  )

  nbasis <- length(knots) + 2

  # mean1  <- c(0, 0, 0, 0, 1, 2, 1, 0,-1, 2, 2,-1, 0, 0.5, 1, 0.5, 0, 0, 0, 0)
  # mean2  <- rev(mean1)
  mean1 <- stats::rnorm(nbasis)
  mean2 <- rev(mean1)
  coef <- rep(0.3, nbasis)

  x1 <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
  x2 <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
  xbind <- rbind(x1$data, x2$data)
  x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
  y <- c(rep(0, n1), rep(1, n2))
  # nu=abs(rowSums(x$x.miss,na.rm=TRUE)/100) # POSSIBLE NU THAT IS DATA RELATED (BAD RESULTS IN 5 ITERATIONS)
  # y=rbinom(n1+n2,1,nu)
  if (unif.miss == TRUE) {
    x <- data.agg.miss(xbind)
  } else {
    print("Beta")
    x <- data.agg.miss.betadecay(xbind, 1)
  }

  data <- data.frame(y = y)
  data["x"] <- x$x.miss
  data <- addgrid(data, grid)
  data
}


#' Data generator for the 2-dimensional partially observed functional data case
#'
#' @param N .
#' @param px .
#' @param py .
#' @param Rsq .
#' @param n_missing .
#' @param min_distance_x .
#' @param min_distance_y .
#'
#' @return .
#'
#' @noRd
data_generator_po2d <- function(N = 100, px = 20, py = 20, Rsq = 0.95, n_missing = 1, min_distance_x = 9, min_distance_y = 9) {
  Data_H <- function(x_observations, y_observations, epsilon_1 = 0.2, epsilon_2 = 0.2, epsilon_data = 0.015) {
    x_b <- length(x_observations)
    y_b <- length(y_observations)

    DATA_T <- DATA_N <- matrix(nrow = x_b, ncol = y_b)

    a1 <- stats::rnorm(1, 0, epsilon_1)
    a2 <- stats::rnorm(1, 0, epsilon_2)

    for (i in 1:x_b) {
      for (j in 1:y_b) {
        DATA_T[i, j] <- a1 * cos(2 * pi * x_observations[i]) + a2 * cos(2 * pi * y_observations[j]) + 1
        DATA_N[i, j] <- DATA_T[i, j] + stats::rnorm(1, 0, epsilon_data)
      }
    }

    DATA <- data.frame(DATA_T[, 1])

    DATA[["DATA_T"]] <- DATA_T
    DATA[["DATA_N"]] <- DATA_N

    DATA <- DATA[, -1]


    # a is the stochastic part of the surface
    # this is needed for the function `response_int_H` !!!!!!!!!!!
    res <- list(
      DATA = DATA,
      a = c(a1, a2),
      epsilon_data = epsilon_data
    )
  }
  Stochastic_Data_H <- function(x, y, a1, a2, epsilon_data = 0.015) {
    # 'e' stands for epsilon

    x_b <- length(x)
    y_b <- length(y)

    DATA_T <- DATA_N <- matrix(nrow = x_b, ncol = y_b)

    for (i in 1:x_b) {
      for (j in 1:y_b) {
        DATA_T[i, j] <- a1 * cos(2 * pi * x[i]) + a2 * cos(2 * pi * y[j]) + 1
        DATA_N[i, j] <- DATA_T[i, j] + stats::rnorm(1, 0, epsilon_data)
      }
    }

    DATA <- data.frame(DATA_T[, 1])

    DATA[["DATA_T"]] <- DATA_T
    DATA[["DATA_N"]] <- DATA_N

    DATA <- DATA[, -1]
    DATA
  }
  Beta_H_saddle <- function(x_observations, y_observations) {
    aux_beta <- matrix(nrow = length(x_observations), ncol = length(y_observations))

    for (i in 1:length(x_observations)) {
      for (h in 1:length(y_observations)) {
        aux_beta[i, h] <- (4 * x_observations[i] - 2)^3 - 3 * (4 * x_observations[i] - 2) * (4 * y_observations[h] - 2)^2
      }
    }
    aux_beta / 10
  }
  Beta_H_exp <- function(x_observations, y_observations) {
    aux_beta <- matrix(nrow = length(x_observations), ncol = length(y_observations))

    for (i in 1:length(x_observations)) {
      for (h in 1:length(y_observations)) {
        aux_beta[i, h] <- 5 * exp(-8 * ((x_observations[i] - 0.75)^2 + (y_observations[h] - 0.75)^2)) + 5 * exp(-8 * ((x_observations[i] - 0.1)^2 + (y_observations[h] - 0.1)^2))
      }
    }
    aux_beta / 10
  }
  response_int_H <- function(f_X, a1, a2, epsilon_data, f_Beta, x_observations, y_observations, sub_response = 50, N = FALSE) {
    n_y <- m_y <- 2 * sub_response

    W_delta_y <- array(dim = (n_y + 1) * (m_y + 1))

    h <- (x_observations[length(x_observations)] - x_observations[1]) / n_y # if y_b and y_a are functions of x_observations
    HX <- (y_observations[length(y_observations)] - y_observations[1]) / m_y # this line should go inside the next for loop (the int_i lopp or outer loop)


    x <- seq(x_observations[1], x_observations[length(x_observations)], h)
    y <- seq(y_observations[1], y_observations[length(y_observations)], HX)


    simp_w_y <- rep(1, n_y + 1)
    even_y <- seq(2, n_y + 1 - 1, 2)
    odd_y <- seq(3, n_y + 1 - 1, 2)

    simp_w_y[even_y] <- 2
    simp_w_y[odd_y] <- 4

    Sim_w_x_y <- (h / 3) * simp_w_y

    W_x_y <- (HX / 3) * Sim_w_x_y
    W_x_even_y <- 2 * W_x_y
    W_x_odd_y <- 4 * W_x_y

    for (aux in 1:(m_y + 1)) {
      # print(c("aux = ", aux))

      if (aux == 1 || aux == (m_y + 1)) {
        W_delta_y[((n_y + 1) * (aux - 1) + 1):((n_y + 1) * aux)] <- W_x_y
      } else {
        if (aux %% 2 == 0) {
          W_delta_y[((n_y + 1) * (aux - 1) + 1):((n_y + 1) * aux)] <- W_x_even_y
        } else {
          W_delta_y[((n_y + 1) * (aux - 1) + 1):((n_y + 1) * aux)] <- W_x_odd_y
        }
      }
    }

    X_eval <- f_X(x, y, a1, a2, epsilon_data)

    if (N) {
      X_eval <- X_eval$DATA_N
    } else {
      X_eval <- X_eval$DATA_T
    }


    Beta_eval <- f_Beta(x, y)

    y_int <- as.double(t(ks::vec(X_eval)) %*% diag(W_delta_y) %*% ks::vec(Beta_eval))

    list(
      nu = y_int,
      X_eval = X_eval
    )
  }


  x <- seq(from = 0, to = 1, length.out = px)
  y <- seq(from = 0, to = 1, length.out = py)

  X <- lapply(1:N, function(a) {
    Data_H(x, y)
  })


  X_true <- lapply(seq_along(X), function(i) X[[i]]$DATA$DATA_T)
  X_real <- lapply(seq_along(X), function(i) X[[i]]$DATA$DATA_N)


  nu <- sapply(1:N, function(a) {
    response_int_H(Stochastic_Data_H,
                   X[[a]]$a[[1]],
                   X[[a]]$a[[2]],
                   X[[a]]$epsilon_data,
                   f_Beta = Beta_H_saddle,
                   x_observations = x, y_observations = y
    )[["nu"]]
  })

  X_miss = add_miss(X_real, n_missing = n_missing, min_distance_x = min_distance_x, min_distance_y = min_distance_y)

  var_e <- (1 / Rsq - 1) * stats::var(nu) # (1-Rsq)*var(nu[ind,])

  response <- nu + stats::rnorm(N, sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL #rnorm(nu[ind,])

  list(
    y = response,
    nu = nu,
    X_true = X_true,
    X_real = X_real,
    X_miss = X_miss[["X_miss"]],
    miss_points = X_miss[["miss_points"]],
    missing_points = X_miss[["missing_points"]]
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
add_miss <- function(X, n_missing = 1, min_distance_x = 9, min_distance_y = 9) {
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
    x_missing <- sort(sample(1:(x_b - min_distance_x), n_missing))
    y_missing <- sort(sample(1:(y_b - min_distance_y), n_missing))

    x_pos <- x_missing:(x_missing + min_distance_x - 1)
    y_pos <- y_missing:(y_missing + min_distance_y - 1)

    missing_points[[j]] <- expand.grid(x_pos, y_pos)

    for (add_miss in 1:nrow(missing_points[[j]])) {
      X[[j]][missing_points[[j]][add_miss, 1], missing_points[[j]][add_miss, 2]] <-
        NA
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
