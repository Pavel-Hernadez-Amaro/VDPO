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
bspline <-function(X., XL., XR., NDX., BDEG.){
  dx <- (XR. - XL.)/NDX.
  knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
  B <- splines::spline.des(knots, X., BDEG.+1, 0*X.)$design
  res <- list(B = B, knots = knots)
  res
}
