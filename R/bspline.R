#' B-spline generator
#'
#' @param X. Points where you are going to evaluate the function.
#' @param XL. Min x (- epsilon)
#' @param XR. Max x (+ epsilon)
#' @param NDX. Number of basis - BDEG.
#' @param BDEG. Degree of the polynomial.
#'
#' @return A \code{list} where the first element is the design matrix of the
#' B-spline and the second element is a list with the different knots.
#'
#' @noRd
bspline <-function(x, xl, xr, nseg, bdeg){
  dx <- (xr - xl)/nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- splines::spline.des(knots, x, bdeg + 1, 0 * x)$design
  list(B = B, knots = knots)
}
