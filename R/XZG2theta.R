#' XZG2theta function
#'
#' @param X Fixed part of the mixed model.
#' @param Z Random part of the mixed model.
#' @param G Variance-covariance matrix.
#' @param TMatrix Matrix of transformation from the multivariate model to the mixed model.
#' @param y Response variable.
#' @param family Family of the distribution.
#' @param offset .
#'
#' @return .
#'
#' @noRd
XZG2theta <- function(X, Z, G, TMatrix, y, family = stats::gaussian(), offset = NULL) {
  if (dim(X)[1] != dim(Z)[1]) {
    stop("'X' and 'Z'must have same numbers of rows", call. = FALSE)
  }

  # if ( dim(Z)[2]!=length(G[[1]]) || dim(Z)[2]!=length(G[[2]]) ) {
  #   stop("The number of columns of 'Z' must be equal to the length of 'G'")
  # }

  w <- as.vector(rep(1, dim(X)[1]))

  fit <- sop.fit(
    X = X, Z = Z, G = G,
    y = y, family = family,
    control = list(trace = FALSE), offset = offset
  )

  if (dim(fit$Vp)[1] - dim(TMatrix)[1] == 0) {
    theta_aux <- c(fit$b.fixed, fit$b.random)
  } else {
    ## nacho:: de la manera en la que esto está programado,
    ## siempre agregar primero la variable no funcional y
    ## luego la funcional

    aux <- dim(fit$Vp)[1] - dim(TMatrix)[1]
    theta_aux <- c(fit$b.fixed[-(1:aux)], fit$b.random) # ESTE [1:4] FUE AGREGADO PARA QUE AL AGREGAR OTRAS VARIABLES LINEALES SOLO COGIERA LA PARTE FUNCIONAL
  }

  theta <- TMatrix %*% theta_aux

  if (dim(fit$Vp)[1] == dim(TMatrix)[1]) {
    covar_theta <- TMatrix %*% fit$Vp %*% t(TMatrix)
    std_error_theta <- sqrt(diag(TMatrix %*% fit$Vp %*% t(TMatrix)))
    std_error_non_functional <- NULL
    p_values <- NULL
  } else {
    covar_theta <- TMatrix %*% fit$Vp[-(1:aux), -(1:aux)] %*% t(TMatrix)
    std_error_theta <- sqrt(diag(TMatrix %*% fit$Vp[-(1:aux), -(1:aux)] %*% t(TMatrix)))
    std_error_non_functional <- sqrt(diag(fit$Vp[1:aux, 1:aux]))
    WALD <- (fit$b.fixed[(1:aux)] / std_error_non_functional)
    p_values <- 2 * stats::pnorm(abs(WALD), lower.tail = FALSE)
  }

  list(
    fit                      = fit,
    theta                    = theta,
    std_error_theta          = std_error_theta,
    std_error_non_functional = std_error_non_functional,
    covar_theta              = covar_theta,
    p_values                 = p_values
  )
}

#' One dimensional version of XZG2theta.
#'
#' @param X Fixed part of the mixed model.
#' @param Z Random part of the mixed model.
#' @param G Variance-covariance matrix.
#' @param TMatrix Matrix of transformation from the multivariate model to the mixed model.
#' @param y Response variable.
#' @param family Family of the distrbution.
#'
#' @return .
#'
#' @noRd
XZG2theta_1d <- function(X, Z, G, TMatrix, y, family = stats::gaussian()) {
  fit <- sop.fit(
    X = X, Z = Z, G = G,
    y = y, family = family, weights = rep(1, dim(X)[1]),
    control = list(trace = FALSE)
  )

  theta_aux <- c(fit$b.fixed, fit$b.random)
  theta <- TMatrix %*% theta_aux

  list(fit = fit, theta = theta)
}
