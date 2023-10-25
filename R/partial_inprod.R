#' Partial inner product
#'
#' @param n_intervals Number of intervals where we are going to integrate.
#' This should be an even number.
#' @param knots1 First set of nodes.
#' @param knots2 Seconds set of nodes.
#' @param bdeg Degree of the basis. This should be a two dimensional vector.
#' @param spline_domain Domain where the spline is defined.
#' @param rng Integration limits. This should be a two dimensional vector.
#'
#' @return Matrix with the value of the integral.
#' @export
partial_inprod <- function(n_intervals, knots1, knots2, bdeg, spline_domain, rng) {
  if (n_intervals %% 2 != 0) {
    stop("The number of intervals should be even", call. = FALSE)
  }

  if (length(rng) != 2) {
    stop("Length should be a vector with two elements", call. = FALSE)
  }

  if (length(bdeg) != 2) {
    stop("bdeg should be a vector with two elements", call. = FALSE)
  }

  width <- (rng[[2]] - rng[[1]]) / n_intervals


  fx <- splines::spline.des(knots1, rng, bdeg[[1]] + 1, 0 * rng)$design
  ft <- splines::spline.des(knots2, rng, bdeg[[2]] + 1, 0 * rng)$design

  fT <- spline_domain
  fBeta <- kronecker(ft, fT)

  XI0 <- matrix(crossprod(ft, fBeta), nrow = ncol(ft), ncol = ncol(fBeta))
  XI1 <- 0
  XI2 <- 0

  for (i in 1:(n_intervals-1)) {
    x <- rng[[1]] + i * width

    fx <- splines::spline.des(knots1, x, bdeg[[1]] + 1, 0 * x)$design
    ft <- splines::spline.des(knots2, x, bdeg[[2]] + 1, 0 * x)$design

    fBeta <- kronecker(ft, fT)

    FX <- matrix(crossprod(fT, fBeta), nrow = ncol(fT), ncol = ncol(fBeta))

    if (i %% 2 == 0) {
      XI2 <- XI2 + FX
    } else {
      XI1 <- XI1 + FX
    }
  }

  XI <- width * (XI0 + 2 * XI2 + 4 * XI1) / 3

  XI
}
