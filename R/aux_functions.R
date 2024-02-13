#' @references All credits to the \href{https://cran.r-project.org/web/packages/sommer/index.html}{sommer} package authors.
Rten2 <- function(X1, X2) {
  one.1 <- matrix(1, 1, ncol(X1))
  one.2 <- matrix(1, 1, ncol(X2))
  kronecker(X1, one.2) * kronecker(one.1, X2)
}

#' @references All credits to the \href{https://cran.r-project.org/web/packages/MASS/}{MASS} package authors.
ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  #
  # based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
  #
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) {
    stop("'X' must be a numeric or complex matrix")
  }
  if (!is.matrix(X)) X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) {
    Xsvd$v %*% (1 / Xsvd$d * t(Xsvd$u))
  } else if (!any(Positive)) {
    array(0, dim(X)[2L:1L])
  } else {
    Xsvd$v[, Positive, drop = FALSE] %*% ((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
  }
}

#' @references This is a modified version of sop.fit to improve numerical stability.
#' All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
construct.matrices <- function(X, Z, z, w) {
  XtW. <- t(X * w)
  XtX. <- XtW. %*% X
  XtZ. <- XtW. %*% Z
  ZtX. <- t(XtZ.)
  ZtW. <- t(Z * w)
  ZtZ. <- ZtW. %*% Z
  Xty. <- XtW. %*% z
  Zty. <- ZtW. %*% z
  yty. <- sum((z^2) * w)
  ZtXtZ <- rbind(XtZ., ZtZ.)
  u <- c(Xty., Zty.)
  res <- list(
    XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ.,
    Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ,
    u = u
  )
}

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
construct.block <- function(A1, A2, A3, A4) {
  block <- rbind(cbind(A1, A2), cbind(A3, A4))
  return(block)
}

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
sop.control <- function(maxit = 200, epsilon = 1e-06, trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0) {
    stop("value of 'epsilon' must be > 0")
  }
  if (!is.numeric(maxit) || maxit <= 0) {
    stop("maximum number of iterations must be > 0")
  }
  list(maxit = maxit, epsilon = epsilon, trace = trace)
}

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
sop.fit <- function(y, X, Z, weights = NULL, G = NULL, vcstart = NULL,
                    etastart = NULL, mustart = NULL, offset = NULL, family = stats::gaussian(),
                    control = SOP::sop.control()) {
  deviance <- function(C, G, w, sigma2, ssr, edf) {
    log_det_C <- determinant(C)$modulus
    log_det_G <- determinant(G)$modulus
    deviance <- log_det_C + log_det_G + sum(log(sigma2 *
      1 / w)) + ssr / sigma2 + edf
    deviance
  }
  control <- do.call("sop.control", control)
  if (missing(X)) {
    stop("Missing argument: 'X' must be provided")
  }
  if (missing(y)) {
    stop("Missing argument: 'y' must be provided")
  }
  if (missing(Z)) {
    stop("Missing argument: 'Z' must be provided")
  }
  if (!is.null(vcstart)) {
    if (length(vcstart) != (length(G) + 1)) {
      stop("The length of 'vcstart' should be equal to the length of 'G' + 1)")
    }
  }
  trace <- control$trace
  ynames <- if (is.matrix(y)) {
    rownames(y)
  } else {
    names(y)
  }
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- NCOL(X)
  EMPTY <- nvars == 0
  nc <- ncol(X)
  ncomp <- length(G)
  nq <- ncol(Z)
  if (!is.null(vcstart)) {
    la <- vcstart
  } else {
    la <- rep(1, len = ncomp + 1)
  }
  devold <- 1e+10
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
  }
  prior.weights <- weights
  if (is.null(offset)) {
    offset <- rep.int(0, nobs)
  }
  if (is.null(G)) {
    stop("Missing argument: 'G' must be provided")
  }
  Xnames <- colnames(X)
  na.ind <- is.na(y)
  y.tmp <- y
  y[na.ind] <- 1
  weights <- weights * (!na.ind)
  X <- as.matrix(X)
  start <- NULL
  variance <- family$variance
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) {
    stop("illegal `family' argument")
  }
  valideta <- family$valideta
  if (is.null(valideta)) {
    valideta <- function(eta) TRUE
  }
  validmu <- family$validmu
  if (is.null(validmu)) {
    validmu <- function(mu) TRUE
  }
  if (is.null(mustart)) {
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (NCOL(y) > 1) {
    stop("y must be univariate")
  }
  eta <- if (!is.null(etastart)) {
    etastart
  } else if (!is.null(start)) {
    if (length(start) != nvars) {
      stop(gettextf(
        "Length of start should equal %d and correspond to initial coefs.",
        nvars
      ))
    } else {
      coefold <- start
      offset + as.vector(if (nvars == 1) {
        c(X, Z) * start
      } else {
        c(X, Z) %*% start
      })
    }
  } else {
    family$linkfun(mustart)
  }
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta))) {
    stop("Can't find valid starting values: please specify some")
  }
  for (i in 1:control$maxit) {
    deriv <- family$mu.eta(eta)
    z <- (eta - offset) + (y - mu) / deriv
    w <- as.vector(deriv^2 / family$variance(mu))
    w <- w * weights
    z[!weights] <- 0
    mat <- construct.matrices(X, Z, z, w) # GLAM = FALSE is now the only option
    if (trace) {
      start1 <- proc.time()[3]
    }
    for (it in 1:control$maxit) {
      Ginv <- 0
      for (ii in 1:ncomp) {
        Ginv <- Ginv + 1 / la[ii + 1] * G[[ii]]
      }
      GG <- 1 / Ginv
      V <- construct.block(
        mat$XtX., mat$XtZ., mat$ZtX.,
        mat$ZtZ.
      )
      D <- diag(c(rep(0, nc), Ginv))
      H <- (1 / la[1]) * V + D
      Hinv <- try(solve(H))
      if (inherits(Hinv, "try-error")) {
        Hinv <- ginv(H)
      }
      b <- (1 / la[1]) * Hinv %*% mat$u
      b.fixed <- b[1:nc]
      b.random <- b[-(1:nc)]
      aux <- GG - diag(Hinv[-(1:nc), -(1:nc)])
      ed.sum <- 0
      ied <- taus <- NULL
      for (i.fun in 1:ncomp) {
        d1.d <- (1 / la[i.fun + 1]) * G[[i.fun]]
        ed1 <- sum(aux * d1.d)
        ed1 <- ifelse(ed1 <= 1e-10, 1e-10, ed1)
        tau1 <- sum(b.random^2 * G[[i.fun]]) / ed1
        tau1 <- ifelse(tau1 <= 1e-10, 1e-10, tau1)
        taus <- c(taus, tau1)
        ied <- c(ied, ed1)
      }
      ssr <- mat$yty. - t(c(b.fixed, b.random)) %*% (2 *
        mat$u - V %*% b)
      dev <- deviance(
        H, diag(GG), w[w != 0], la[1], ssr,
        sum(b.random^2 * Ginv)
      )[1]
      if (family$family == "gaussian" | family$family ==
        "Gamma" | family$family == "quasipoisson") {
        sig2 <- as.numeric((ssr / (length(y[weights !=
          0]) - sum(ied) - nc)))
      } else {
        sig2 <- 1
      }
      lanew <- c(sig2, taus)
      dla <- abs(devold - dev)
      if (trace) {
        cat(sprintf("%1$3d %2$10.6f", it, dev))
        cat(sprintf("%8.3f", ied), "\n")
      }
      if (dla < control$epsilon) {
        break
      }
      la <- lanew

      if (la[1] < 1e-9) { ####### this if is new
        la[1] <- 1e-9
      }

      devold <- dev
    }
    if (trace) {
      end1 <- proc.time()[3]
      cat("Timings:\nSOP", (end1 - start1), "seconds\n")
    }
    eta.old <- eta
    eta <- X %*% b.fixed + Z %*% b.random + offset
    mu <- linkinv(eta)
    tol <- sum((eta - eta.old)^2) / sum(eta^2)
    if (tol < control$epsilon | (family$family == "gaussian" &
      family$link == "identity")) {
      break
    }
  }
  end <- proc.time()[3]
  mu.eta <- family$mu.eta
  mu.eta.val <- mu.eta(eta)
  linear.predictor <- eta
  mu <- linkinv(eta)
  names(mu) <- ynames
  names(linear.predictor) <- ynames
  residuals <- family$dev.resids(y, mu, weights)
  s <- attr(residuals, "sign")
  if (is.null(s)) {
    s <- sign(y - mu)
  }
  residuals <- sqrt(pmax(residuals, 0)) * s
  names(residuals) <- ynames
  residuals[na.ind] <- NA
  names(b.fixed) <- Xnames
  names(b.random) <- colnames(Z)
  names(la) <- c("ssr", names(G))
  names(ied) <- c(names(G))
  dev.residuals <- family$dev.resids(y, mu, weights)
  dev.residuals[na.ind] <- NA
  deviance <- sum(dev.residuals, na.rm = TRUE)
  null.deviance <- stats::glm(y ~ offset(offset),
    family = family,
    weights = prior.weights
  )$deviance
  out <- list(
    tol.ol = tol, it.ol = i, tol.il = dla, it.in = it,
    vc = la, edf = ied
  )
  fit <- list()
  fit$b.fixed <- b.fixed
  fit$b.random <- b.random
  fit$fitted.values <- mu
  fit$linear.predictor <- linear.predictor
  fit$residuals <- residuals
  fit$X <- X
  fit$Z <- Z
  fit$G <- G
  fit$y <- y.tmp
  fit$weights <- weights
  fit$family <- family
  fit$out <- out
  fit$deviance <- deviance
  fit$null.deviance <- null.deviance
  fit$Vp <- Hinv
  class(fit) <- "sop"
  invisible(fit)
}

#' Data generator function
#'
#' @param N Number of subjects.
#' @param J Number of maximum observations per subject.
#' @param nsims Number of simulations per the simulation study.
#' @param aligned If the data is aligned or not.
#' @param multivariate If TRUE, the data is generated with 2 variables.
#' @param Rsq .
#' @param use_x If the data is generated with x.
#' @param use_f If the data is generated with f.
#'
#' @return Example data.
#'
#' @noRd
dg <- function(N = 100, J = 100, nsims = 1, Rsq = 0.95, aligned = TRUE, multivariate = FALSE, use_x = FALSE, use_f = FALSE) {
  for (iter in 1:nsims) {
    set.seed(42 + iter)

    # Generating the domain for all subject with a minimum of 10 observations
    M <- round(stats::runif(N, 10, J), digits = 0)

    if (max(M) > J) {
      M[which(M > J)] <- J
    }

    if (min(M) <= 10) {
      M[which(M <= 10)] <- 10
    }

    maxM <- max(M)
    t <- 1:maxM

    # We can sort the data without loss of generality
    M <- sort(M)

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

      Beta[i, 1:(M[i]), 1] <- ((10 * t[1:(M[i])] / M[i]) - 5) / 10
      Beta[i, 1:(M[i]), 2] <- ((1 - (2 * M[i] / maxM)) * (5 - 40 * ((t[1:(M[i])] / M[i]) - 0.5)^2)) / 10

      if (multivariate) {
        nu[i] <- sum(X_s[i, ] * Beta[i, , 1], na.rm = 1) / (M[i]) + sum(Y_s[i, ] * Beta[i, , 2], na.rm = 1) / (M[i])
      } else {
        nu[i] <- sum(X_s[i, ] * Beta[i, , 1], na.rm = 1) / (M[i]) # NOT NOISY
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

#' add grid for ffpo
#'
#' This function should be only used when the \code{bidimensional_grid}
#' parameter of \code{ffpo} is \code{FALSE}.
#'
#' @param df .
#' @param grid .
#'
#' @return Dataframe with grid
#' @export
#'
#' @noRd
addgrid <- function(df, grid) {
  N <- nrow(df)
  l <- length(grid)

  newgrid <- suppressWarnings(matrix(grid, nrow = N))
  newgrid[(l + 1):length(newgrid)] <- NA

  df["grid"] <- newgrid
  df
}

Data_H <- function(x_observations, y_observations, epsilon_1 = 0.2, epsilon_2 = 0.2, epsilon_data = 0.015) {
  x_b <- length(x_observations)
  y_b <- length(y_observations)

  DATA_T <- DATA_N <- matrix(nrow = x_b, ncol = y_b)

  a1 <- rnorm(1, 0, epsilon_1)
  a2 <- rnorm(1, 0, epsilon_2)

  for (i in 1:x_b) {
    for (j in 1:y_b) {
      DATA_T[i, j] <- a1 * cos(2 * pi * x_observations[i]) + a2 * cos(2 * pi * y_observations[j]) + 1
      DATA_N[i, j] <- DATA_T[i, j] + rnorm(1, 0, epsilon_data)
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
      DATA_N[i, j] <- DATA_T[i, j] + rnorm(1, 0, epsilon_data)
    }
  }

  DATA <- data.frame(DATA_T[, 1])

  DATA[["DATA_T"]] <- DATA_T
  DATA[["DATA_N"]] <- DATA_N

  DATA <- DATA[, -1]
  DATA
}

x <- y <- seq(from = 0, to = 1, length.out = 20)
Data_H(x, y) -> res
plotly::plot_ly(z = res$DATA_T, type = "surface")
plotly::plot_ly(z = res$DATA_N, type = "surface")


#' Title
#'
#' @param x_observations
#' @param y_observations
#'
#' @return
#'
#' @examples
Beta_fun <- function(x_observations, y_observations) {
  x_b <- length(x_observations)
  y_b <- length(y_observations)

  Beta_x <- -1 * (5 - 40 * ((x_observations / x_b) - 0.5)^2) / 50
  Beta_y <- -1 * (5 - 40 * ((y_observations / y_b) - 0.5)^2) / 50

  res <- matrix(0, nrow = length(x_observations), ncol = length(y_observations))
  res <- kronecker(t(Beta_y), Beta_x)
  res
}
Beta_fun(x, y)

Beta_H_saddle <- function(x_observations, y_observations) {
  aux_beta <- matrix(nrow = length(x_observations), ncol = length(y_observations))

  for (i in 1:length(x_observations)) {
    for (h in 1:length(y_observations)) {
      aux_beta[i, h] <- (4 * x_observations[i] - 2)^3 - 3 * (4 * x_observations[i] - 2) * (4 * y_observations[h] - 2)^2
    }
  }
  aux_beta / 10
}

# plot_ly(z = Beta_H_saddle(x, y), type = "surface")


Beta_H_exp <- function(x_observations, y_observations) {
  aux_beta <- matrix(nrow = length(x_observations), ncol = length(y_observations))

  for (i in 1:length(x_observations)) {
    for (h in 1:length(y_observations)) {
      aux_beta[i, h] <- 5 * exp(-8 * ((x_observations[i] - 0.75)^2 + (y_observations[h] - 0.75)^2)) + 5 * exp(-8 * ((x_observations[i] - 0.1)^2 + (y_observations[h] - 0.1)^2))
    }
  }
  aux_beta / 10
}

# plot_ly(z = Beta_H_exp(x, y), type = "surface")

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

# response_int_H(Data_H, Beta_H_saddle, x, y, N = FALSE)


