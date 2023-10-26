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

  x <- seq(rng[[1]], rng[[2]], width)

  fx <- splines::spline.des(knots1, x, bdeg[[1]] + 1, 0 * x)$design
  ft <- splines::spline.des(knots2, x, bdeg[[2]] + 1, 0 * x)$design

  fT <- spline_domain
  fBeta <- kronecker(ft, fT)

  w <- seq_along(x)
  aux_1 <- ifelse(w %% 2 == 0,
                  2,
                  4)
  aux_1[1] <- 1
  aux_1[length(aux_1)] <- 1
  W <- diag(aux_1)

  W <- width * W / 3

  t(fx) %*% W %*% fBeta

}
