#' Title
#'
#' @param X .
#' @param M .
#' @param nbasis .
#' @param bdeg .
#' @param sub .
#' @param lim .
#'
#' @return .
#' @export
#'
Data2B_simpson <- function(X, M, nbasis = c(30,30,30), bdeg = c(3,3,3), sub = 25, lim=NULL){
  ## SETTING SOME MATRICES AND PARAMETERS

  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1

  N <- nrow(X)

  if (length(M)!=N) {
    stop("length of 'M' must be equal to 'N'")
  }

  c1 <- nbasis[1]
  c2 <- nbasis[2]
  c3 <- nbasis[3]

  # error=K=NULL

  error <- NULL
  K     <- NULL

  rng <- matrix(0, ncol = 2, nrow = N)

  # L_Phi_aux=L_X_aux=L_Phi=L_X=L_y=L_theta=list()
  # vector(mode = "list", length = N)

  L_Phi_aux <- list()
  L_X_aux   <- list()
  L_Phi     <- list()
  L_X       <- list()
  L_y       <- list()
  L_theta   <- list()

  A <- matrix(0,nrow = N, ncol = N*c1)

  for (i in 1:N) {
    if (length(X[i,]) - length(which(is.na(X[i,]))) != M[i]) {
      stop("Incorrect numbers of NAs in column",i, " of 'X'")
    }

    ############### HERE WE CREATE THE BASIS FOR THE DATA

    rng[i,] <- c(1, M[i])

    # XL=rng[i,1]-0.001
    # XR=rng[i,2]+0.001

    XL <- rng[i,1] - 1e-6
    XR <- rng[i,2] + 1e-6

    c <- c1 - bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    # nam_X <- paste("B", i, sep = "_")
    # L_X[[i]]=assign(nam_X, bspline(1:M[i], XL, XR, c, bdeg[1]))

    L_X[[i]] <- bspline(1:M[i], XL, XR, c, bdeg[1])

    # matplot(L_X[[i]]$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT

    ######### ESTIMATING THE DATA COEFFICIENTS (MATRIX A)

    aux <- L_X[[i]]$B
    aux_2 <- B2XZG_1d(aux,2,c1)
    aux_3 <- XZG2theta_1d(X = aux_2$X, Z = aux_2$Z, G = aux_2$G, T = aux_2$T, y = X[i,1:M[i]])

    A[i,((c1*(i-1))+1):(i*c1)] <- aux_3$theta

    L_theta[[i]] <- aux_3$theta

    # nam_y <- paste("y_h", i, sep = "_")
    # L_y[[i]]=assign(nam_y, L_X[[i]]$B%*%L_theta[[i]])

    L_y[[i]] <- L_X[[i]]$B %*% L_theta[[i]]
    error[i] <- mean(abs((X[i,1:M[i]]) - L_y[[i]])) # THE ERROR OF SMOOTHING THE DATA

    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)

    c_t <- c2 - bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    # nam_Phi <- paste("Phi", i, sep = "_")
    # L_Phi[[i]]=assign(nam_Phi, bspline(1:M[i], XL, XR, c_t, bdeg[2]))

    L_Phi[[i]] <- bspline(1:M[i], XL, XR, c_t, bdeg[2])
  }

  ####### HERE WE CREATE THE MARGINAL BASIS FOR THE T VARIABLE In B(t,T)

  xlim_T <- c(min(M), max(M))

  # XL_T=xlim_T[1]-0.001
  # XR_T=xlim_T[2]+0.001

  XL_T <- xlim_T[1] - 1e-06
  XR_T <- xlim_T[2] + 1e-06

  c_T <- c3 - bdeg[3] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

  # if (is.null(lim)) {
  #
  #   B_T <- bspline(M, XL_T, XR_T, c_T, bdeg[3])
  #
  # }else{
  #
  #   B_T <- bspline(M, lim[1]-0.001, lim[2]+0.001, c_T, bdeg[3])
  # }

  B_T <- ifelse(is.null(lim),
                bspline(M, XL_T, XR_T, c_T, bdeg[3]),
                bspline(M, lim[1] - 1e-06, lim[2] + 1e-06, c_T, bdeg[3])
         )


  # matplot(B_T$B,type="l") ## THIS PLOT WILL HELP TO SEE IF THE B-SPLINES BASIS ARE CORRECT

  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT
  # IS NOT MECESSARY TO DO THIS FOR THE T MARGINAL BASIS
  for (i in 1:N) {

    # DATA BASIS

    breaks <- L_X[[i]]$knots

    dife   <- diff(breaks)[1]
    breaks <- c(breaks[1] - dife, breaks, breaks[length(breaks)] + dife)
    breaks <- c(breaks[1] - dife, breaks, breaks[length(breaks)] + dife)
    n      <- length(breaks)

    # nam_X_aux <- paste("B_aux", i, sep = "_")
    # L_X_aux[[i]]=assign(nam_X_aux, create.bspline.basis(breaks=breaks,norder=bdeg[1]+1,dropin=c(1:5,(n-2):(n+2))))

    L_X_aux[[i]] <- fda::create.bspline.basis(breaks = breaks, norder = bdeg[1] + 1, dropin = c(1:5,(n-2):(n+2)))

    # t MARGINAL BASIS

    breaks <- L_Phi[[i]]$knots

    dife   <- diff(breaks)[1]
    breaks <- c(breaks[1] - dife,breaks, breaks[length(breaks)] + dife)
    breaks <- c(breaks[1] - dife,breaks, breaks[length(breaks)] + dife)
    n      <- length(breaks)

    # nam_Phi_aux <- paste("Phi_aux", i, sep = "_")
    # L_Phi_aux[[i]]=assign(nam_Phi_aux, create.bspline.basis(breaks=breaks,norder=bdeg[2]+1,dropin=c(1:5,(n-2):(n+2))))

    L_Phi_aux[[i]] <- fda::create.bspline.basis(breaks = breaks, norder = bdeg[2] + 1, dropin = c(1:5,(n-2):(n+2)))

  }

  # PERFORMING THE INNER PRODUCT


  # need to rewrite this for statement
  for (i in 1:N) {
    PROD <- Simpson(L_X_aux[[i]], L_Phi_aux[[i]], B_T$B[i, ], rng = c(1,M[i]), sub = sub) / M[i]
    K    <- rbind(K,PROD)
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
