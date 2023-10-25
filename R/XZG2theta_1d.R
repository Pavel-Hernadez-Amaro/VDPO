#' One dimensional version of XZG2theta.
#'
#' @param X Fixed part of the mixed model.
#' @param Z Random part of the mixed model.
#' @param G Variance-covariance matrix.
#' @param T Matrix of transformation from the multivariate model to the mixed model.
#' @param y Response variable.
#' @param family Family of the distrbution.
#'
#' @return .
#' @export
XZG2theta_1d <- function(X, Z, G, TMatrix, y, family = stats::gaussian()){
  fit <- sop.fit(
    X = X, Z = Z, G = G,
    y = y, family = family, weights = rep(1,dim(X)[1]),
    control = list(trace = FALSE)
  )

  theta_aux <- c(fit$b.fixed,fit$b.random)
  theta     <- TMatrix %*% theta_aux

  list(fit = fit, theta = theta)
}
