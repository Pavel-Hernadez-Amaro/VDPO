#' data2B simpson for classification problems
#'
#' @param X .
#' @param M .
#' @param grid .
#' @param nbasis .
#' @param bdeg .
#'
#' @return .
#' @export
data2B_simpson_classification <- function(X, M, grid, nbasis=c(80,80), bdeg=c(3,3)){

  # X is the matrix of Data # dim(X) == N x max(M_b)
  # M is the matrix with columns M_a and M_b of numbers of observations dim(M) == N x 2

  # Both X and M should be matrices

  if (!is.matrix(X)) {
    stop("argument 'X' should be a matrix", call. = FALSE)
  }

  if (!is.matrix(M)) {
    stop("argument 'M' should be a matrix", call. = FALSE)
  }

  sub <- 1000

  N <- nrow(X)

  c1 <- nbasis[1]
  c2 <- nbasis[2]

  M_a <- M[,1]
  M_b <- M[,2]
  error <- NULL
  K <- NULL

  rng <- matrix(0,ncol = 2,nrow = N)

  L_Phi     <- vector(mode = "list", length = N)
  L_Phi_aux <- vector(mode = "list", length = N)
  L_X       <- vector(mode = "list", length = N)
  L_y       <- vector(mode = "list", length = N)
  L_theta   <- vector(mode = "list", length = N)


  A <- matrix(0,nrow = N, ncol = N*c1)


  for (i in 1:N) {
    if (length(X[i,]) - length(which(is.na(X[i,])))!=M_b[i]-M_a[i]+1) {
      stop("incorrect numbers of NAs in column ",i, " of 'X'", call. = FALSE)
    }

    # Here we create the basis for the data

    rng[i,] <- c(grid[M_a[i]], grid[M_b[i]]) #c(range(grid)[1],range(grid)[2])

    XL <- rng[i,1] - 1e-06
    XR <- rng[i,2] + 1e-06

    c <- c1 - bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    L_X[[i]] <- bspline(grid[M_a[i]:M_b[i]], XL, XR, c, bdeg[1])


    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)

    aux <- L_X[[i]]$B

    aux_2 <- B2XZG_1d(aux,2,c1)

    aux_3 <- XZG2theta_1d(X = aux_2$X, Z = aux_2$Z, G = aux_2$G, TMatrix = aux_2$T, y=X[i,M_a[i]:M_b[i]] )

    A[i,((c1*(i-1))+1):(i*c1)] <- aux_3$theta

    L_theta[[i]] <- aux_3$theta
    L_y[[i]] <- L_X[[i]]$B%*%L_theta[[i]]

    error[i] <- mean(abs((X[i,M_a[i]:M_b[i]]) - L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA

  }

  # HERE WE CREATE THE BASIS FOR THE t VARIABLE In B(t)

  rng_t <- range(grid) # THE KNOTS FOR THE FUNCTIONAL COEFFICIENTS HAVE TO BE MAXIMUM

  XL_t <- rng_t[1] - 1e-06
  XR_t <- rng_t[2] + 1e-06

  c_t <- c2 - bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

  Phi <- bspline(grid, XL_t, XR_t, c_t, bdeg[2])



  # PERFORMING THE INNER PRODUCT

  res <- matrix(nrow = N, ncol = c2)

  for (i in 1:N) {


    Phi_short <- Phi$B[M_a[i]:M_b[i],]

    aux_del <- which(colSums(Phi_short) == 0)

    if (length(aux_del)!=0)
      Phi_short <- as.matrix(Phi_short[,-aux_del])

    min_knot <- max(which(Phi$knots<=range(grid[M_a[i]:M_b[i]],na.rm=1)[1])) - 3
    max_knot <- min(which(Phi$knots>=range(grid[M_a[i]:M_b[i]],na.rm=1)[2])) + 3

    if (min_knot<1)
      min_knot <- 1

    if (max_knot > length(Phi$knots))
      max_knot <- length(Phi$knots)

    new_knots <- Phi$knots[min_knot:max_knot]

    B_X_a   <- splines::spline.des(L_X[[i]]$knots, rng[i,1], bdeg[1]+1, 0*rng[i,1])$design
    B_Phi_a <- splines::spline.des(new_knots, rng[i,1], bdeg[2]+1, 0*rng[i,1])$design

    B_X_b   <- splines::spline.des(L_X[[i]]$knots, rng[i,2], bdeg[1]+1, 0*rng[i,2])$design
    B_Phi_b <- splines::spline.des(new_knots, rng[i,2], bdeg[2]+1, 0*rng[i,2])$design

    n <- 2*sub # Intervals for the Simpson method

    width <- (rng[i,2] - rng[i,1])/n

    XIa <- t(B_X_a) %*% B_Phi_a
    XIb <- t(B_X_b) %*% B_Phi_b
    XI1 <- 0
    XI2 <- 0

    for (j in 1:(n-1)) {
      x <- rng[i,1] + j*width
      fx1 <- splines::spline.des(L_X[[i]]$knots, x, bdeg[1]+1, 0*x)$design
      fx2 <- splines::spline.des(new_knots, x, bdeg[2]+1, 0*x)$design
      Fx <- t(fx1) %*% fx2

      if (j%%2 == 0) {
        XI2 <- XI2 + Fx
      } else {
        XI1 <- XI1 + Fx
      }
    }

    XI <- matrix( 0,nrow = c1, ncol = c2)

    if (length(aux_del) != 0) {
      XI[,-aux_del] <- width*(XIa + XIb + 2*XI2 + 4*XI1)/3
    } else {
      XI <- width*(XIa + XIb + 2*XI2 + 4*XI1)/3
    }

    res[i,] <- t(L_theta[[i]]) %*% XI

    # This is the same as filling a full-size inner product
    # matrix with zeros and the multiplying it by the matrix A
  }

  list(
    B=res,
    A=A,
    K=K,
    x_h=L_y,
    error=error,
    B_X=L_X,
    theta_x=L_theta,
    Phi=Phi
  )
}

