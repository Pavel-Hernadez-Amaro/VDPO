#' Partial inner product
#'
#' @param n_intervals number of intervals where we are going to integrate.
#' This should be an even number.
#' @param knots1 first set of nodes.
#' @param knots2 seconds set of nodes.
#' @param bdeg degree of the basis. This should be a two dimensional vector.
#' @param spline_domain domain where the spline is defined.
#' @param rng integration limits. This should be a two dimensional vector.
#'
#' @return matrix with the value of the integral.
#'
#' @noRd
partial_inprod <- function(n_intervals, knots1, knots2, bdeg, spline_domain, rng) {
  if (n_intervals %% 2 != 0) {
    stop("the 'n_intervals' parameter should be an even number", call. = FALSE)
  }

  if (length(rng) != 2) {
    stop("'rng' should be a vector with two elements", call. = FALSE)
  }

  if (length(bdeg) != 2) {
    stop("'bdeg' should be a vector with two elements", call. = FALSE)
  }

  width <- (rng[[2]] - rng[[1]]) / n_intervals

  x <- seq(rng[[1]], rng[[2]], width)

  fx <- splines::spline.des(knots1, x, bdeg[[1]] + 1, 0 * x)$design
  ft <- splines::spline.des(knots2, x, bdeg[[2]] + 1, 0 * x)$design

  fT <- spline_domain
  fBeta <- kronecker(ft, fT)

  w <- seq_along(x)
  aux_1 <- ifelse(w %% 2 == 0, 4, 2)
  aux_1[1] <- 1
  aux_1[length(aux_1)] <- 1
  W <- diag(aux_1)

  W <- width * W / 3

  t(fx) %*% W %*% fBeta
}
