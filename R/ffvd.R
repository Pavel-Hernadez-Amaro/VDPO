#' Defining variable domain functional data terms in vd_fit formulae
#'
#' Auxiliary function used to define \code{ffvd} terms within \code{vd_fit} model formulae.
#' This term represents a functional predictor where each function is observed over a domain of varying length.
#' The formulation is \eqn{\frac{1}{T_i} \int _1^{T_i} X_i(t)\beta(t,T_i)dt}, where \eqn{X_i(t)} is a functional covariate of length \eqn{T_i}, and \eqn{\beta(t,T_i)} is an unknown bivariate functional coefficient.
#' The functional basis used to model this term is the B-spline basis.
#'
#' @param X variable domain functional covariate \code{matrix}.
#' @param grid observation points of the variable domain functional covariate.
#' If not provided, it will be `1:ncol(X)`.
#' @param nbasis number of bspline basis to be used.
#' @param bdeg degree of the bspline basis used.
#'
#' @return the function is interpreted in the formula of a \code{VDPO} model.
#' \code{list} containing the following elements:
#' - An item named `B` design matrix.
#' - An item named `X_hat` smoothed functional covariate.
#' - An item named `L_Phi` and \code{B_T} 1-dimensional marginal B-spline basis used for the functional coefficient.
#' - An item named `M` matrix object indicating the observed domain of the data.
#' - An item named `nbasis` number of basis used.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' data <- data_generator_vd(beta_index = 1, use_x = FALSE, use_f = FALSE)
#' X <- data$X_se
#'
#' # Specifying a custom grid
#' custom_grid <- seq(0, 1, length.out = ncol(X))
#' ffvd_term_custom_grid <- ffvd(X, grid = custom_grid)
#'
#' # Customizing the number of basis functions
#' ffvd_term_custom_basis <- ffvd(X, nbasis = c(10, 10, 10))
#'
#' # Customizing both basis functions and degrees
#' ffvd_term_custom <- ffvd(X, nbasis = c(10, 10, 10), bdeg = c(3, 3, 3))
#'
#' @export
ffvd <- function(X, grid, nbasis = c(30, 50, 30), bdeg = c(3, 3, 3)) {

  if (length(dim(X)) != 2) {
    stop("'X' should be 2-dimensional", call. = FALSE)
  }

  if (missing(grid)) {
    grid <- 1:ncol(X)
  }

  sub <- 500
  pord <- c(2, 2)

  # X is the matrix of Data # dim(X) == N x max(M)
  # grid is the vector of observation points
  # M is the vector of numbers of observations dim(M) == N x 2

  M <- t(apply(X, 1, function(x) range(which(!is.na(x)))))

  if (any(M[, 1] >= M[, 2])) {
    stop("no curve can have a negative number of observations", call. = FALSE)
  }

  N <- nrow(X)

  X_hat <- matrix(nrow = N, ncol = ncol(X))

  c1 <- nbasis[1]
  c2 <- nbasis[2]
  c3 <- nbasis[3]

  K <- NULL
  rng <- cbind(grid[M[,1]],grid[M[,2]])

  L_X <- vector(mode = "list", length = N)

  A <- matrix(0, nrow = N, ncol = N * c1)

  for (i in 1:N) {
    ############### HERE WE CREATE THE BASIS FOR THE DATA
    XL <- rng[i, 1] - 1e-6
    XR <- rng[i, 2] + 1e-6

    c <- c1 - bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    L_X[[i]] <- bspline(grid[M[i, 1]:M[i, 2]], XL, XR, c, bdeg[1])

    ######### Estimating the coefficients of the data (matrix A)

    aux <- L_X[[i]]$B
    aux_2 <- B2XZG_1d(aux, pord[1], c1)
    aux_3 <- XZG2theta_1d(X = aux_2$X, Z = aux_2$Z, G = aux_2$G, TMatrix = aux_2$T, y = X[i, M[i, 1]:M[i, 2]])

    X_hat[i, M[i,1]:M[i,2]] <- aux_3$fit$fitted.values

    A[i, ((c1 * (i - 1)) + 1):(i * c1)] <- aux_3$theta

  }

  ####### HERE WE CREATE THE MARGINAL BASIS FOR THE T VARIABLE In B(t,T)

  M_diff <- (M[, 2] - M[, 1] + 1)

  xlim_T <- c(grid[min(M_diff)], grid[max(M_diff)]) ##

  XL_T <- xlim_T[1] - 1e-06
  XR_T <- xlim_T[2] + 1e-06

  c_T <- c3 - bdeg[3] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

  if (all(M[, 1] == 1)) {
    B_T <- bspline(grid[M[, 2]], XL_T, XR_T, c_T, bdeg[3])
  } else {
    B_T <- bspline(grid[(M[, 2] - M[, 1] + 1)], XL_T, XR_T, c_T, bdeg[3])
  }

  ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE in B(t,T)

  c_t <- c2 - bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

  L_X_all <- bspline(grid, min(grid) - 1e-6, max(grid) + 1e-6, c_t, bdeg[2])

  # PERFORMING THE INNER PRODUCT
  for (i in 1:N) {
    PROD <- partial_inprod(
      n_intervals   = sub,
      knots1        = L_X[[i]]$knots,
      knots2        = L_X_all$knots,
      bdeg          = bdeg[1:2],
      spline_domain = B_T$B[i, , drop = FALSE],
      rng           = c(grid[M[i, 1]], grid[M[i, 2]])
    )
    PROD <- PROD / grid[(M[i, 2] - M[i, 1] + 1)]

    K <- rbind(K, PROD)
  }

  B <- A %*% K

  res <- list(
    B      = B,
    X_hat  = X_hat,
    L_Phi  = L_X_all,
    B_T    = B_T,
    M      = M,
    nbasis = nbasis
  )

  class(res) <- "ffvd"

  res
}



#' @noRd
B2XZG <- function(B_all, deglist) {
  nffvd <- length(deglist)

  pord <- c(2, 2)

  Tnlist <- vector(mode = "list", length = nffvd)
  Tslist <- vector(mode = "list", length = nffvd)

  t1list <- vector(mode = "list", length = nffvd)
  t2list <- vector(mode = "list", length = nffvd)

  G <- vector(mode = "list", length = 2 * nffvd)

  for (i in seq_along(deglist)) {
    c1 <- deglist[[i]][1]
    c2 <- deglist[[i]][2]

    D_1 <- diff(diag(c1), differences = pord[1])
    D_2 <- diff(diag(c2), differences = pord[2])

    P1.svd <- svd(crossprod(D_1))
    P2.svd <- svd(crossprod(D_2))

    U_1s <- P1.svd$u[, 1:(c1 - pord[1])] # eigenvectors
    U_1n <- P1.svd$u[, -(1:(c1 - pord[1]))]
    d1 <- P1.svd$d[1:(c1 - pord[1])] # eigenvalues

    U_2s <- P2.svd$u[, 1:(c2 - pord[2])] # eigenvectors
    U_2n <- P2.svd$u[, -(1:(c2 - pord[2]))]
    d2 <- P2.svd$d[1:(c2 - pord[2])] # eigenvalues

    Tnlist[[i]] <- kronecker(U_1n, U_2n) # this is T_n

    AUX_1 <- kronecker(U_1n, U_2s)
    AUX_2 <- kronecker(U_1s, U_2n)
    AUX_3 <- kronecker(U_1s, U_2s)

    Tslist[[i]] <- cbind(AUX_1, AUX_2, AUX_3) # this is T_s

    d_1s <- diag(P1.svd$d)[1:(c1 - pord[1]), 1:(c1 - pord[1])]
    d_2s <- diag(P2.svd$d)[1:(c2 - pord[2]), 1:(c2 - pord[2])]

    T_1 <- kronecker(diag(pord[1]), d_2s)
    T_2 <-
      matrix(0,
        nrow = pord[2] * (c1 - pord[1]),
        ncol = pord[2] * (c1 - pord[1])
      )
    T_3 <- kronecker(diag(c1 - pord[1]), d_2s)

    T_21 <-
      cbind(T_1, matrix(
        0,
        nrow = dim(T_1)[1],
        ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_1)[2]
      ))
    T_22 <-
      cbind(
        matrix(0, nrow = dim(T_2)[1], ncol = dim(T_1)[2]),
        T_2,
        matrix(
          0,
          nrow = dim(T_2)[1],
          ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_1)[2] - dim(T_2)[2]
        )
      )
    T_23 <-
      cbind(matrix(
        0,
        nrow = ((c2 - pord[2]) * (c1 - pord[1])),
        ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_3)[2]
      ), T_3)

    H_1 <- matrix(0,
      nrow = pord[1] * (c2 - pord[2]),
      ncol = pord[1] * (c2 - pord[2])
    )
    H_2 <- kronecker(d_1s, diag(pord[2]))
    H_3 <- kronecker(d_1s, diag(c2 - pord[2]))

    H_11 <-
      cbind(H_1, matrix(
        0,
        nrow = dim(H_1)[1],
        ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_1)[2]
      ))
    H_12 <-
      cbind(
        matrix(0, nrow = dim(H_2)[1], ncol = dim(H_1)[2]),
        H_2,
        matrix(
          0,
          nrow = dim(H_2)[1],
          ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_1)[2] - dim(H_2)[2]
        )
      )
    H_13 <-
      cbind(matrix(
        0,
        nrow = ((c2 - pord[2]) * (c1 - pord[1])),
        ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_3)[2]
      ), H_3)

    L_2 <- rbind(T_21, T_22, T_23)
    L_1 <- rbind(H_11, H_12, H_13)

    t1list[[i]] <- diag(L_2) # no mistake, diag(L_2) is t1
    t2list[[i]] <- diag(L_1)
  }

  T_n <- Matrix::bdiag(Tnlist)
  T_s <- Matrix::bdiag(Tslist)

  TMatrix <- cbind(T_n, T_s)

  Z <- B_all %*% T_s
  X <- B_all %*% T_n

  it <- 1
  for (i in seq_along(deglist)) {
    G[[2 * i - 1]] <- rep(0, ncol(Z))
    G[[2 * i]] <- rep(0, ncol(Z))

    G[[2 * i - 1]][it:((it + length(t1list[[i]])) - 1)] <- t1list[[i]]
    G[[2 * i]][it:((it + length(t2list[[i]])) - 1)] <- t2list[[i]]

    it <- it + length(t1list[[i]])
  }

  TMatrix <- cbind(T_n, T_s)

  list(
    X_ffvd  = as.matrix(X),
    Z_ffvd  = as.matrix(Z),
    G_ffvd  = G,
    TMatrix = as.matrix(TMatrix)
  )
}
