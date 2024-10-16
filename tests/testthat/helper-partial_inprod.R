partial_inprod_arguments_generator <- function() {
  sub <- 500
  pord <- c(2, 2)

  VDPO_example_vd <- data_generator_vd()

  X <- VDPO_example_vd$X_se
  M <- rng <- t(apply(X, 1, function(x) range(which(!is.na(x)))))
  N <- nrow(X)

  c1 <- c2 <- c3 <- 10
  bdeg <- c(3, 3, 3)

  L_Phi <- vector(mode = "list", length = N)
  L_X <- vector(mode = "list", length = N)

  A <- matrix(0, nrow = N, ncol = N * c1)
  for (i in 1:N) {
    XL <- rng[i, 1] - 1e-6
    XR <- rng[i, 2] + 1e-6

    c <- c1 - bdeg[1]

    L_X[[i]] <- bspline(M[i, 1]:M[i, 2], XL, XR, c, bdeg[1])

    aux <- L_X[[i]]$B
    aux_2 <- B2XZG_1d(aux, pord[1], c1)
    aux_3 <- XZG2theta_1d(X = aux_2$X, Z = aux_2$Z, G = aux_2$G, TMatrix = aux_2$T, y = X[i, M[i, 1]:M[i, 2]])

    A[i, ((c1 * (i - 1)) + 1):(i * c1)] <- aux_3$theta

    c_t <- c2 - bdeg[2]

    L_Phi[[i]] <- bspline(M[i, 1]:M[i, 2], XL, XR, c_t, bdeg[2])
  }

  M_diff <- (M[, 2] - M[, 1] + 1)

  xlim_T <- c(min(M_diff), max(M_diff))

  XL_T <- xlim_T[1] - 1e-06
  XR_T <- xlim_T[2] + 1e-06
  c_T <- c3 - bdeg[3]

  if (all(M[, 1] == 1)) {
    B_T <- bspline(M[, 2], XL_T, XR_T, c_T, bdeg[3])
  } else {
    B_T <- bspline((M[, 2] - M[, 1] + 1), XL_T, XR_T, c_T, bdeg[3])
  }

  # case i=1
  list(
    n_intervals   = sub,
    knots1        = L_X[[1]]$knots,
    knots2        = L_Phi[[1]]$knots,
    bdeg          = bdeg[1:2],
    spline_domain = B_T$B[1, , drop = FALSE],
    rng           = c(M[1, 1], M[1, 2])
  )
}
