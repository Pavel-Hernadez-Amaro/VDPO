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

#' B-spline generator
#'
#' @param X. Points where you are going to evaluate the function.
#' @param XL. Min x (- epsilon).
#' @param XR. Max x (+ epsilon).
#' @param NDX. Number of basis - BDEG.
#' @param BDEG. Degree of the polynomial.
#'
#' @return A \code{list} where the first element is the design matrix of the
#' B-spline and the second element is a list with the different knots.
#'
#' @noRd
bspline <- function(x, xl, xr, nseg, bdeg) {
  dx <- (xr - xl) / nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- splines::spline.des(knots, x, bdeg + 1, 0 * x)$design
  list(B = B, knots = knots)
}

#' This is a modified version of sop.fit to improve numerical stability.
#'
#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
#'
#' @noRd
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
construct.block <- utils::getFromNamespace("construct.block", "SOP")

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
sop.control <- utils::getFromNamespace("sop.control", "SOP")

#' Slightly modified version of \code{sop.fit}. We have asked the authors of the package for permission to make this change.
#'
#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
#'
#' @noRd
sop.fit <- function(y, X, Z, weights = NULL, G = NULL, vcstart = NULL,
                    etastart = NULL, mustart = NULL, offset = NULL, family = stats::gaussian(),
                    control = sop.control()) {
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

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
RH <- utils::getFromNamespace("RH", "SOP")

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
construct.1D.pspline <- utils::getFromNamespace("construct.1D.pspline", "SOP")

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
construct.2D.pspline <- utils::getFromNamespace("construct.2D.pspline", "SOP")

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
construct.3D.pspline <- utils::getFromNamespace("construct.3D.pspline", "SOP")

#' @references All credits to the \href{https://cran.r-project.org/web/packages/SOP/index.html}{SOP} package authors.
construct.capital.lambda <- utils::getFromNamespace("construct.capital.lambda", "SOP")

#' Bidimensional functional data generator with random errors.
#'
#' @param x_observations Observations of the x axis.
#' @param y_observations Observations of the y axis.
#' @param epsilon_1,epsilon_2 Standard deviation of the stochastic components.
#' @param epsilon_data Standard deviation of the noise.
#'
#' @return Bidimensional functional data.
#'
#' @noRd
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

#' Bidimensional functional data generator with fixed errors.
#'
#' @param x Observations of the x axis.
#' @param y Observations of the x axis.
#' @param a1,a2 Fixed parameters of the stochastic components.
#' @param epsilon_data Standard deviation of the noise.
#'
#' @return Bidimensional functional data.
#'
#' @noRd
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

#' Bidimensional functional coefficient using 2-dimensional b-splines.
#'
#' @param x_observations Observations of the x axis.
#' @param y_observations Observations of the y axis.
#'
#' @return Bidimensional functional coefficient.
#'
#' @noRd
Beta_fun <- function(x_observations, y_observations) {
  x_b <- length(x_observations)
  y_b <- length(y_observations)

  Beta_x <- -1 * (5 - 40 * ((x_observations / x_b) - 0.5)^2) / 50
  Beta_y <- -1 * (5 - 40 * ((y_observations / y_b) - 0.5)^2) / 50

  res <- matrix(0, nrow = length(x_observations), ncol = length(y_observations))
  res <- kronecker(t(Beta_y), Beta_x)
  res
}

#' Bidimensional functional coefficient with saddle shape.
#'
#' @param x_observations Observations of the x axis.
#' @param y_observations Observations of the y axis.
#'
#' @return Bidimensional functional coefficient.
#'
#' @noRd
Beta_H_saddle <- function(x_observations, y_observations) {
  aux_beta <- matrix(nrow = length(x_observations), ncol = length(y_observations))

  for (i in 1:length(x_observations)) {
    for (h in 1:length(y_observations)) {
      aux_beta[i, h] <- (4 * x_observations[i] - 2)^3 - 3 * (4 * x_observations[i] - 2) * (4 * y_observations[h] - 2)^2
    }
  }
  aux_beta / 10
}

#' Bidimensional functional coefficient with exponential shape.
#'
#' @param x_observations Observations of the x axis.
#' @param y_observations Observations of the y axis.
#'
#' @return Bidimensional functional coefficient.
#'
#' @noRd
Beta_H_exp <- function(x_observations, y_observations) {
  aux_beta <- matrix(nrow = length(x_observations), ncol = length(y_observations))

  for (i in 1:length(x_observations)) {
    for (h in 1:length(y_observations)) {
      aux_beta[i, h] <- 5 * exp(-8 * ((x_observations[i] - 0.75)^2 + (y_observations[h] - 0.75)^2)) + 5 * exp(-8 * ((x_observations[i] - 0.1)^2 + (y_observations[h] - 0.1)^2))
    }
  }
  aux_beta / 10
}

#' Response generator for functional regression models.
#'
#' @param f_X,f_Beta Functions to be integrated.
#' @param a1,a2 Fixed parameters of the stochastic components.
#' @param epsilon_data Standard deviation of the noise for the covariates.
#' @param x_observations Observations of the x axis.
#' @param y_observations Observations of the y axis.
#' @param sub_response Number of intervals for the Simpson integration method.
#' @param noise Boolean indicating the presence of noise. Defaults to \code{FALSE}.
#'
#' @return Response variable for the functional regression model.
#'
#' @noRd
response_int_H <- function(f_X, a1, a2, epsilon_data, f_Beta, x_observations, y_observations, sub_response = 50, noise = FALSE) {
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

  if (noise) {
    X_eval <- X_eval$DATA_N
  } else {
    X_eval <- X_eval$DATA_T
  }


  Beta_eval <- f_Beta(x, y)

  y_int <- as.double(t(vec(X_eval)) %*% diag(W_delta_y) %*% vec(Beta_eval))

  list(
    nu = y_int,
    X_eval = X_eval
  )
}

#' From matrix to vector
#'
#' @references All credits to the \href{https://cran.r-project.org/web/packages/ks/index.html}{ks} package authors.
#'
#' @noRd
vec <- function(x, byrow = FALSE) {
  if (is.vector(x))
    return(x)
  if (byrow)
    x <- t(x)
  d <- ncol(x)
  vecx <- vector()
  for (j in 1:d) vecx <- c(vecx, x[, j])
  return(vecx)
}

#' From vector to matrix
#'
#' @references All credits to the \href{https://cran.r-project.org/web/packages/ks/index.html}{ks} package authors.
#'
#' @noRd
invvec <- function (x, ncol, nrow, byrow = FALSE) {
  if (length(x) == 1)
    return(x)
  d <- sqrt(length(x))
  if (missing(ncol) | missing(nrow)) {
    ncol <- d
    nrow <- d
    if (round(d) != d)
      stop("Need to specify nrow and ncol for non-square matrices")
  }
  invvecx <- matrix(0, nrow = nrow, ncol = ncol)
  if (byrow)
    for (j in 1:nrow) invvecx[j, ] <- x[c(1:ncol) + (j -
                                                       1) * ncol]
  else for (j in 1:ncol) invvecx[, j] <- x[c(1:nrow) + (j -
                                                          1) * nrow]
  return(invvecx)
}

#' @noRd
dummy <- function() {
  SOP::f
  utils::getFromNamespace
}

#' Helper function to add zeros to the side of a vector
#'
#' @noRd
add_zeros_to_side <- function(vector, final_length, side = c("right", "left")) {
  side <- match.arg(side)

  zeros <- rep(0, final_length - length(vector))

  if (side == "right") {
    return(c(vector, zeros))
  } else if (side == "left"){
    return(c(zeros, vector))
  }
}

#' Helper function to convert a list with vectors and matrices to a data.frame
#' where the matrices are represented as a single columns
#'
#' @noRd
list_to_df <- function(data, response) {
  res <- data.frame(data[[response]])
  names(res) <- toString(response)

  for (variable in setdiff(names(data), response)) {
    res[[variable]] <- data[[variable]]
  }

  res
}



