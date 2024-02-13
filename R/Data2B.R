#' Data to B matrix
#'
#' @param X .
#' @param M Cantidad de elementos que hay en las filas de X que no es NA.
#' Se puede suprimir.
#' @param nbasis .
#' @param bdeg .
#' @param sub .
#' @param lim .
#'
#' @return .
#' @export
Data2B <- function(X, M, nbasis = c(30, 30, 30), bdeg = c(3, 3, 3), sub = 500, lim = NULL) {
  ## SETTING SOME MATRICES AND PARAMETERS

  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1

  N <- nrow(X)

  if (length(M) != N) {
    stop("length of 'M' must be equal to the number of rows of X.")
  }

  c1 <- nbasis[1]
  c2 <- nbasis[2]
  c3 <- nbasis[3]

  error <- NULL
  K <- NULL

  rng <- matrix(0, ncol = 2, nrow = N)

  L_Phi <- vector(mode = "list", length = N)
  L_X <- vector(mode = "list", length = N)
  L_y <- vector(mode = "list", length = N)
  L_theta <- vector(mode = "list", length = N)

  A <- matrix(0, nrow = N, ncol = N * c1)

  for (i in 1:N) {
    if (length(X[i, ]) - length(which(is.na(X[i, ]))) != M[i]) {
      stop(paste0("Incorrect numbers of NAs in column", i, " of 'X'"),
        call. = FALSE
      )
    }

    ############### HERE WE CREATE THE BASIS FOR THE DATA

    rng[i, ] <- c(1, M[i])

    XL <- rng[i, 1] - 1e-6
    XR <- rng[i, 2] + 1e-6

    c <- c1 - bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1


    L_X[[i]] <- bspline(1:M[i], XL, XR, c, bdeg[1])

    ######### Estimating the coefficients of the data (matrix A)

    aux <- L_X[[i]]$B
    aux_2 <- B2XZG_1d(aux, 2, c1)
    aux_3 <- XZG2theta_1d(X = aux_2$X, Z = aux_2$Z, G = aux_2$G, TMatrix = aux_2$T, y = X[i, 1:M[i]])

    A[i, ((c1 * (i - 1)) + 1):(i * c1)] <- aux_3$theta

    L_theta[[i]] <- aux_3$theta

    L_y[[i]] <- L_X[[i]]$B %*% L_theta[[i]]
    error[i] <- mean(abs((X[i, 1:M[i]]) - L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA

    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)

    c_t <- c2 - bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    L_Phi[[i]] <- bspline(1:M[i], XL, XR, c_t, bdeg[2])
  }

  ####### HERE WE CREATE THE MARGINAL BASIS FOR THE T VARIABLE In B(t,T)

  xlim_T <- c(min(M), max(M))

  XL_T <- xlim_T[1] - 1e-06
  XR_T <- xlim_T[2] + 1e-06

  c_T <- c3 - bdeg[3] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

  if (is.null(lim)) {
    B_T <- bspline(M, XL_T, XR_T, c_T, bdeg[3])
  } else {
    B_T <- bspline(M, lim[1] - 0.001, lim[2] + 0.001, c_T, bdeg[3])
  }
  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT
  # IS NOT MECESSARY TO DO THIS FOR THE T MARGINAL BASIS

  # PERFORMING THE INNER PRODUCT


  # need to rewrite this for statement
  for (i in 1:N) {
    PROD <- partial_inprod(
      n_intervals = sub,
      knots1 = L_X[[i]]$knots,
      knots2 = L_Phi[[i]]$knots,
      bdeg = bdeg[1:2],
      spline_domain = B_T$B[i, , drop = FALSE],
      rng = c(1, M[i])
    )
    PROD <- PROD / M[i]

    K <- rbind(K, PROD)
  }

  res <- A %*% K

  list(
    B       = res,
    A       = A,
    K       = K,
    x_h     = L_y,
    error   = error,
    B_X     = L_X,
    theta_x = L_theta,
    B_T     = B_T,
    B_Phi   = L_Phi
  )
}
