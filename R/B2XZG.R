#' Multidimensional B2XZG function
#'
#' @param B Design matrix.
#' @param pord Order of the matrix of differences. We are going to use a
#' bidimensional penalization.
#' @param c Polynomial order.
#'
#' @return .
#' @export
#'
B2XZG <-function (B, pord = c(2,2),c = c(10,10)) {

  c1 <- c[1]
  c2 <- c[2]
  c1c2 <-  ncol(B)

  if (c1c2 != c1 * c2) {
    stop("c1 * c2 must me equal to the number of colums of B", call. = FALSE)
  }

  D_1 <- diff(diag(c1), differences = pord[1])
  D_2 <- diff(diag(c2), differences = pord[2])

  P1.svd <- svd(crossprod(D_1))
  P2.svd <- svd(crossprod(D_2))

  U_1s <- (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
  U_1n <- ((P1.svd$u)[,-(1:(c1-pord[1]))])
  d1   <- (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues

  U_2s <- (P2.svd$u)[,1:(c2-pord[2])] # eigenvectors
  U_2n <- ((P2.svd$u)[,-(1:(c2-pord[2]))])
  d2   <- (P2.svd$d)[1:(c2-pord[2])]  # eigenvalues


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

  T <- cbind(T_n,T_s)

  ####

  list(
    X    = X,
    Z    = Z,
    G    = G,
    T    = T,
    d1   = d1,
    d2   = d2,
    D_1  = D_1,
    D_2  = D_2,
    U_1n = U_1n,
    U_1s = U_1s,
    U_2n = U_2n,
    U_2s = U_2s,
    T_n  = T_n,
    T_s  = T_s,
    t_1  = t_1,
    t_2  = t_2
  )
}
