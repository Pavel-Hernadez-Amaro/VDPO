## NO ALINEAR A LA IZDA LAS CURVAS. LINEAS A CAMBIAR CON ##

ffvd <- function(X, nbasis = c(30, 30, 30), bdeg = c(3, 3, 3)) {
  sub <- 500
  pord <- c(2,2)
  ## SETTING SOME MATRICES AND PARAMETERS

  # X is the matrix of Data # dim(X) == N x max(M)
  # M is the vector of numbers of observations dim(M) == N x 1

  M <- t(apply(X, 1, function(x) range(which(!is.na(x)))))

  if (any(M[, 1] >= M[, 2]))
    stop("Ninguna curva puede tener un número de observaciones negativa", call. = FALSE)

  N <- nrow(X)

  c1 <- nbasis[1]
  c2 <- nbasis[2]
  c3 <- nbasis[3]

  K <- NULL
  rng <- M

  L_Phi <- vector(mode = "list", length = N)
  L_X   <- vector(mode = "list", length = N)

  A <- matrix(0, nrow = N, ncol = N * c1)

  for (i in 1:N) {

    ############### HERE WE CREATE THE BASIS FOR THE DATA
    XL <- rng[i, 1] - 1e-6
    XR <- rng[i, 2] + 1e-6

    c <- c1 - bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1


    L_X[[i]] <- bspline(M[i, 1]:M[i, 2], XL, XR, c, bdeg[1])

    ######### Estimating the coefficients of the data (matrix A)

    aux   <- L_X[[i]]$B
    aux_2 <- B2XZG_1d(aux, pord[1], c1)
    aux_3 <- XZG2theta_1d(X = aux_2$X, Z = aux_2$Z, G = aux_2$G, TMatrix = aux_2$T, y = X[i,M[i, 1]:M[i, 2]])

    A[i,((c1 * (i - 1)) + 1):(i * c1)] <- aux_3$theta

    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)

    c_t <- c2 - bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    L_Phi[[i]] <- bspline(M[i, 1]:M[i, 2], XL, XR, c_t, bdeg[2])
  }

  ####### HERE WE CREATE THE MARGINAL BASIS FOR THE T VARIABLE In B(t,T)

  M_diff= (M[,2]-M[,1]+1)

  xlim_T <- c(min(M_diff), max(M_diff)) ##

  XL_T <- xlim_T[1] - 1e-06
  XR_T <- xlim_T[2] + 1e-06

  c_T <- c3 - bdeg[3] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

  if (all(M[,1]==1)) {

    B_T <- bspline(M[,2], XL_T, XR_T, c_T, bdeg[3])

  }else{
    B_T <- bspline((M[,2]-M[,1]+1), XL_T, XR_T, c_T, bdeg[3])
  }
  ################# HERE WE ARE GOING TO TRASNFORM THE B-SPLINES BASIS INTO THE RAMSAY TYPE B-SPLINES BASIS TO PERFORM THE INNER PRODUCT
  # IS NOT MECESSARY TO DO THIS FOR THE T MARGINAL BASIS

  # PERFORMING THE INNER PRODUCT


  # need to rewrite this for statement
  for (i in 1:N) {
    PROD <- partial_inprod(n_intervals   = sub,
                           knots1        = L_X[[i]]$knots,
                           knots2        = L_Phi[[i]]$knots,
                           bdeg          = bdeg[1:2],
                           spline_domain = B_T$B[i, , drop = FALSE],
                           rng           = c(M[i, 1], M[i, 2]))
    PROD <- PROD / (M[i, 2] - M[i, 1] + 1)

    K    <- rbind(K,PROD)
  }

  B <- A %*% K

#### HAY QUE CREAR UNA NUEVA FUNCIÓN A PARTIR DE AQUÍ PARA UN MODELO ADDITIVE

  ##### EMPEZAMOS A TRANSFORMAR EL MODELO MULTIVARIANTE AL MODELO MIXTO (B2XGZ)

  D_1 <- diff(diag(c1), differences = pord[1])
  D_2 <- diff(diag(c2), differences = pord[2])

  P1.svd <- svd(crossprod(D_1))
  P2.svd <- svd(crossprod(D_2))

  U_1s <- P1.svd$u[,1:(c1-pord[1])] # eigenvectors
  U_1n <- P1.svd$u[,-(1:(c1-pord[1]))]
  d1   <- P1.svd$d[1:(c1-pord[1])]  # eigenvalues

  U_2s <- P2.svd$u[,1:(c2-pord[2])] # eigenvectors
  U_2n <- P2.svd$u[,-(1:(c2-pord[2]))]
  d2   <- P2.svd$d[1:(c2-pord[2])]  # eigenvalues

  T_n <- kronecker(U_1n,U_2n)

  AUX_1 <- kronecker(U_1n,U_2s)
  AUX_2 <- kronecker(U_1s,U_2n)
  AUX_3 <- kronecker(U_1s,U_2s)

  T_s <- cbind(AUX_1,AUX_2,AUX_3)

  Z <- B %*% T_s
  X <- B %*% T_n

  ####

  d_1s <- diag(P1.svd$d)[1:(c1-pord[1]),1:(c1-pord[1])]
  d_2s <- diag(P2.svd$d)[1:(c2-pord[2]),1:(c2-pord[2])]

  T_1 <- kronecker(diag(pord[1]),d_2s)
  T_2 <- matrix(0, nrow = pord[2] * (c1-pord[1]), ncol = pord[2] * (c1-pord[1]))
  T_3 <- kronecker(diag(c1-pord[1]),d_2s)

  T_21 <- cbind(T_1, matrix(0, nrow = dim(T_1)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_1)[2]))
  T_22 <- cbind(matrix(0, nrow = dim(T_2)[1], ncol = dim(T_1)[2]), T_2, matrix(0, nrow = dim(T_2)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_1)[2] - dim(T_2)[2]))
  T_23 <- cbind(matrix(0,nrow=((c2-pord[2])*(c1-pord[1])),ncol=(c1*c2-pord[1]*pord[2])-dim(T_3)[2]),T_3)

  H_1 <- matrix(0, nrow = pord[1] * (c2 - pord[2]),ncol = pord[1] * (c2 - pord[2]))
  H_2 <- kronecker(d_1s, diag(pord[2]))
  H_3 <- kronecker(d_1s, diag(c2-pord[2]))

  H_11 <- cbind(H_1, matrix(0, nrow = dim(H_1)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_1)[2]))
  H_12 <- cbind(matrix(0, nrow = dim(H_2)[1],ncol = dim(H_1)[2]), H_2, matrix(0, nrow = dim(H_2)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_1)[2] - dim(H_2)[2]))
  H_13 <- cbind(matrix(0, nrow = ((c2 - pord[2]) * (c1 - pord[1])),ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_3)[2]), H_3)

  L_2 <- rbind(T_21, T_22, T_23)
  L_1 <- rbind(H_11, H_12, H_13)

  t_2 <- diag(L_1)
  t_1 <- diag(L_2)

  G <- list(t_1, t_2)
  names(G) <- c("t_1", "t_2")

  TMatrix <- cbind(T_n,T_s)

  list(
    X_ffvd  = X,
    Z_ffvd  = Z,
    G_ffvd  = G,
    L_Phi   = L_Phi,
    B_T     = B_T,
    M       = M,
    TMatrix = TMatrix
  )


}
