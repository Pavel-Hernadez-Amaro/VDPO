#' One dimensional B2XZG
#'
#' @param B Design matrix.
#' @param pord Order of the difference matrix.
#' @param c Polynomial order.
#'
#' @return .
#'
#' @noRd
B2XZG_1d <- function(B, pord = 2, c = 10) {
  c1 <- c[1]
  D_1 <- diff(diag(c1), differences = pord)
  P1.svd <- svd(crossprod(D_1))

  U_1s <- P1.svd$u[, 1:(c1 - pord[1])] # eigenvectors
  U_1n <- P1.svd$u[, -(1:(c1 - pord[1]))]
  d1 <- P1.svd$d[1:(c1 - pord[1])] # eigenvalues

  T_n <- U_1n
  T_s <- U_1s

  Z <- B %*% T_s
  X <- B %*% T_n

  G <- list(d1)

  T_ <- cbind(T_n, T_s)

  list(
    X    = X,
    Z    = Z,
    G    = G,
    T    = T_,
    d1   = d1,
    D_1  = D_1,
    U_1n = U_1n,
    U_1s = U_1s,
    T_n  = T_n,
    T_s  = T_s
  )
}

#' Special case of the B2XZG_1d function.
#'
#' @param B_all .
#' @param deglist .
#'
#' @return .
#'
#' @noRd
B2XZG_ffpo <- function(B_all, deglist) {
  pord <- 2
  nffpo <- length(deglist)

  Tnlist <- vector(mode = "list", length = nffpo)
  Tslist <- vector(mode = "list", length = nffpo)
  tlist <- vector(mode = "list", length = nffpo)

  G <- vector(mode = "list", length = nffpo)

  for (i in seq_along(deglist)) {
    c1 <- deglist[[i]]
    D_1 <- diff(diag(c1), differences = pord)
    P1.svd <- svd(crossprod(D_1))

    U_1s <- P1.svd$u[, 1:(c1 - pord)] # eigenvectors
    U_1n <- P1.svd$u[, -(1:(c1 - pord))]
    d1 <- P1.svd$d[1:(c1 - pord)] # eigenvalues

    Tslist[[i]] <- U_1s
    Tnlist[[i]] <- U_1n
    tlist[[i]] <- d1
  }

  T_n <- Matrix::bdiag(Tnlist)
  T_s <- Matrix::bdiag(Tslist)

  Z <- B_all %*% T_s
  X <- B_all %*% T_n

  it <- 1
  for (i in seq_along(deglist)) {
    G[[i]] <- rep(0, ncol(Z))
    G[[i]][it:((it + length(tlist[[i]])) - 1)] <- tlist[[i]]

    it <- it + length(tlist[[i]])
  }

  TMatrix <- cbind(T_n, T_s)

  list(
    X = as.matrix(X),
    Z = as.matrix(Z),
    G = G,
    TMatrix = as.matrix(TMatrix)
  )
}
