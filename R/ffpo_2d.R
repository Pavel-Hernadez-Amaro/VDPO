#' Defining partially observed bidimensional functional data terms in VDPO formulae
#'
#' Auxiliary function used to define \code{ffpo_2d} terms within \code{VDPO} model
#' formulae.
#'
#' The difference between
#' miss_points and missing_points is the format in which the data is
#' presented.
#'
#' @param X partially observed bidimensional functional covariate \code{matrix}.
#' @param miss_points,missing_points \code{list} of missing observation points.
#' See 'Details' for more information about the difference in structure between both.
#' @param nbasis number of basis to be used.
#' @param bdeg degree of the basis to be used.
#'
#' @return The function is interpreted in the formula of a \code{VDPO} model.
#' \code{list} containing the following elements:
#' - \code{B_ffpo2d} design matrix.
#' - \code{Phi_ffpo2d} bidimensional B-spline basis used for the functional coefficient.
#' - \code{M_ffpo2d} the \code{missing_points} used as input in the function.
#' - \code{nbasis} number of the basis used.
#'
#' @details
#' \code{miss_points} is a \code{list} of \code{list}s where each inner list corresponds
#' to the observation points in the y-axis and contains the observation points
#' of the missing values for the x-axis. \code{miss_points} acts as a guide for
#' identifying and addressing missing observations in functional data and is used
#' for properly calculating the inner product matrix.
#'
#' \code{missing_points} is a \code{list} where each element is a \code{matrix}
#' containing the missing observations points.
#'
#' @seealso \code{\link{VDPO}}
#'
#' @export
ffpo_2d <- function(X, miss_points, missing_points, nbasis = rep(15, 4), bdeg = rep(3, 4)) {
  for (i in 2:length(X)) {
    if (all(dim(X[[1]]) != dim(X[[i]]))) {
      stop("all matrices inside 'X_miss' should have the same dimensions",
        call. = FALSE
      )
    }
  }

  x_b <- nrow(X[[1]])
  y_b <- ncol(X[[1]])
  N <- length(X)
  x_observations <- seq(from = 0, to = 1, length.out = x_b)
  y_observations <- seq(from = 0, to = 1, length.out = y_b)

  c1 <- nbasis[1]
  c2 <- nbasis[2]
  c1_beta <- nbasis[3]
  c2_beta <- nbasis[4]

  ####### DEFINING SOME EMPTY LISTS, ARRAYS, VECTORS AND MATRICES
  A <- array(0, dim = c(N, N * c1 * c2))
  X_hat <- vector(mode = "list", length = N)

  # plotly::plot_ly(z = X_hat[[1]], type = "surface")

  ####### DEFINING THE PARAMETERS OF THE INNER PRODUCT

  Inner_matrix <- array(dim = c(N * c1 * c2, c1_beta * c2_beta))

  sub <- 500
  n <- 2 * sub # THIS IS THE NUMBER OF INTERVALS FOR THE SIMPSON METHOD
  m <- 2 * sub

  W_delta <- array(dim = (n + 1) * (m + 1))

  simp_w <- rep(1, n + 1)
  even <- seq(2, n + 1 - 1, 2)
  odd <- seq(3, n + 1 - 1, 2)

  simp_w[even] <- 4
  simp_w[odd] <- 2

  h <- (x_observations[x_b] - x_observations[1]) / n # if y_b and y_a are functions of x
  HX <- (y_observations[y_b] - y_observations[1]) / m # this line should go inside the next for loop (the int_i lopp or outer loop)

  Sim_w_x <- (h / 3) * simp_w

  x <- seq(x_observations[1], x_observations[x_b], h)
  y <- seq(y_observations[1], y_observations[y_b], HX)

  W_x <- (HX / 3) * Sim_w_x
  W_x_even <- 4 * W_x
  W_x_odd <- 2 * W_x

  for (aux in 1:(m + 1)) {
    # print(c("aux = ", aux))

    if (aux == 1 || aux == (m + 1)) {
      W_delta[((n + 1) * (aux - 1) + 1):((n + 1) * aux)] <- W_x
    } else {
      if (aux %% 2 == 0) {
        W_delta[((n + 1) * (aux - 1) + 1):((n + 1) * aux)] <- W_x_even
      } else {
        W_delta[((n + 1) * (aux - 1) + 1):((n + 1) * aux)] <- W_x_odd
      }
    }
  }

  W_delta_save <- W_delta


  ####### DEFINING THE B-spline BASIS USED IN THE MODEL (only depend on the observation grid and the number
  # of basis c)

  aux_model_k <-
    bspline(
      x_observations,
      x_observations[1] - 1e-04,
      x_observations[x_b] + 1e-04,
      c1 - bdeg[1],
      bdeg[1]
    ) # x_observations and y_observation have to be the complete domain
  aux_2_model_k <-
    bspline(
      y_observations,
      y_observations[1] - 1e-04,
      y_observations[y_b] + 1e-04,
      c2 - bdeg[2],
      bdeg[2]
    ) # this domain may not by fulfill by any of the surfaces

  knots_x <- aux_model_k$knots
  knots_y <- aux_2_model_k$knots

  fx <- splines::spline.des(knots_x, x, bdeg[1] + 1, 0 * x)$design
  fy <- splines::spline.des(knots_y, y, bdeg[2] + 1, 0 * y)$design

  aux_model <- aux_model_k$B
  aux_2_model <- aux_2_model_k$B

  B_kron_model <- kronecker(aux_2_model, (aux_model))

  aux_model_beta <- bspline(
    x_observations,
    x_observations[1] - 1e-04,
    x_observations[x_b] + 1e-04,
    c1_beta - bdeg[3],
    bdeg[3]
  ) # x_observations and y_observation have to be the complete domain
  aux_2_model_beta <-
    bspline(
      y_observations,
      y_observations[1] - 1e-04,
      y_observations[y_b] + 1e-04,
      c2_beta - bdeg[4],
      bdeg[4]
    ) # this domain may not by fulfill by any of the surfaces

  B_kron_beta <- kronecker(aux_2_model_beta$B, aux_model_beta$B)


  knots_x_beta <- aux_model_beta$knots
  knots_y_beta <- aux_2_model_beta$knots

  fx_beta <- splines::spline.des(knots_x_beta, x, bdeg[3] + 1, 0 * x)$design
  fy_beta <- splines::spline.des(knots_y_beta, y, bdeg[4] + 1, 0 * y)$design

  for (j in 1:N) {
    # print(paste0("j = ", j))

    NA_ind <- NULL
    for (ind_i in 1:dim(X[[j]])[1]) {
      for (ind_j in 1:dim(X[[j]])[2]) {
        if (is.na(X[[j]][ind_i, ind_j])) {
          NA_ind <-
            sort(c(NA_ind, (dim(X[[j]])[1] * (ind_j - 1)) + ind_i))
        }
      }
    }


    I <- diag(dim(B_kron_model)[1])
    I[NA_ind, ] <- 0

    aux_22 <-
      B2XZG_2d(I %*% B_kron_model,
        pord = c(2, 2),
        c = c(c1, c2)
      )

    aux_33 <-
      XZG2theta(
        X = aux_22$X,
        Z = aux_22$Z,
        G = aux_22$G,
        TMatrix = aux_22$TMatrix,
        y = c(X[[j]])
      ) # , weights = w)

    A[j, ((c2 * c1 * (j - 1)) + 1):(j * c1 * c2)] <- aux_33$theta

    X_hat[[j]] <-
      matrix(t(aux_33$theta) %*% t(B_kron_model),
        nrow = x_b,
        ncol = y_b
      )

    for (aux in 1:dim(missing_points[[j]])[1]) {
      X_hat[[j]][missing_points[[j]][aux, 1], missing_points[[j]][aux, 2]] <- NA
    }
    # plotly::plot_ly(z = X_hat[[1]], type = "surface")

    # HERE BEGINS THE DOUBLE INTEGRAL

    W_delta <- W_delta_save

    for (int_i in seq_along(y)) {
      na_x <- NULL

      # print(c("int_i = ", int_i))

      if (length(y) != m + 1) {
        stop("'length(y)' has to be equal to m+1", call. = FALSE)
      }

      prev_obs <- max(which(y[int_i] >= y_observations))

      next_obs <- min(which(y[int_i] <= y_observations))

      na_22 <- unique(append(miss_points[[j]][[prev_obs]], miss_points[[j]][[next_obs]]))

      observed_points <- which(!(seq(1:nrow(X[[j]])) %in% na_22))

      for (aux_o in observed_points) {
        aux_prev <- aux_o - 1
        aux_next <- aux_o + 1

        if ((aux_prev %in% na_22) && (aux_next %in% na_22)) {
          na_22 <- sort(c(na_22, aux_o))
          observed_points <- observed_points[-which(aux_o == observed_points)]
        }
      }

      where_diff <- which(diff(observed_points) != 1)

      if (!(length(where_diff) == 0)) {
        observed_points_final <- NULL

        observed_points_1 <- range(observed_points[1:where_diff[1]])

        if (length(where_diff) > 1) {
          for (aux_diff in seq(length(where_diff) - 1)) {
            observed_points_final <- c(
              observed_points_final,
              range(observed_points[(where_diff[aux_diff] + 1):where_diff[aux_diff + 1]])
            )
          }
        }

        observed_points_last <-
          range(observed_points[(where_diff[length(where_diff)] + 1):length(observed_points)])

        observed_points_group <-
          c(
            observed_points_1,
            observed_points_final,
            observed_points_last
          )
      } else {
        observed_points_group <- range(observed_points)
      }


      ###########

      where_diff_na_22 <- which(diff(na_22) != 1)

      if (!(length(where_diff_na_22) == 0)) {
        na_22_final <- NULL

        na_22_1 <- range(na_22[1:(where_diff_na_22[1])])
        na_22_1[2] <- na_22_1[2] + 1

        if (na_22[1] != 1) {
          na_22_1[1] <- na_22_1[1] - 1
        }

        if (length(where_diff_na_22) > 1) {
          for (aux_diff in seq(length(where_diff_na_22) - 1)) {
            aux <-
              range(na_22[(where_diff_na_22[aux_diff] + 1):where_diff_na_22[aux_diff + 1]])
            aux[1] <- aux[1] - 1
            aux[2] <- aux[2] + 1

            na_22_final <- c(na_22_final, aux)
          }
        }

        na_22_last <-
          c(na_22[length(na_22)] - 1, na_22[length(na_22)])

        if (na_22[length(na_22)] != nrow(X[[j]])) {
          na_22_last[length(na_22_last)] <- na_22_last[length(na_22_last)] + 1
        }

        na_22_group <- c(na_22_1, na_22_final, na_22_last)
      } else {
        if (!(length(na_22) == 0)) {
          na_22_group <- range(na_22)

          if (na_22_group[2] != nrow(X[[j]])) {
            na_22_group[2] <- na_22_group[2] + 1
          }

          if (na_22_group[1] != 1) {
            na_22_group[1] <- na_22_group[1] - 1
          }
        } else {
          na_22_group <- NULL
        }
      }

      if (!is.null(na_22_group)) {
        for (tag_x in 1:(length(na_22_group) / 2)) {
          na_x <-
            c(
              na_x,
              which(x >= x_observations[na_22_group[(2 * tag_x - 1)]] & x <= x_observations[na_22_group[2 * tag_x]])
            )
        }

        W_delta[((n + 1) * (int_i - 1) + 1):((n + 1) * int_i)][na_x] <- 0
      }
    } # for in int_i

    W <- invvec(W_delta, ncol = m + 1, nrow = n + 1)

    # check=10
    # all.equal(as.matrix(W_delta[((n+1)*(check-1)+1):(check*(n+1))]),as.matrix(W[,check]))

    aux_GLAM <- RH(t(Rten2((fy_beta), fy)), RH(t(Rten2((
      fx_beta
    ), fx)), W))
    dim(aux_GLAM) <- c(c1, c1_beta, c2, c2_beta)
    aux_GLAM_apperm <- matrix(aperm(aux_GLAM, c(1, 3, 2, 4)), nrow = c1 * c2)

    # all.equal(aux,t(aux_GLAM))

    # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_apperm))
    # all.equal(Inner_matrix_WO[(c1*c2*(j-1)+1):(c1*c2*j),],(aux_GLAM_apperm))

    Inner_matrix[(c1 * c2 * (j - 1) + 1):(c1 * c2 * j), ] <- aux_GLAM_apperm/nrow(missing_points[[1]])
  }

  B_ffpo2d <- A %*% Inner_matrix

  list(
    B_ffpo2d = B_ffpo2d,
    A = A,
    Inner_matrix = Inner_matrix,
    X_hat=X_hat,
    Phi_ffpo2d = B_kron_beta,
    M_ffpo2d = missing_points,
    nbasis = nbasis
  )
}


#' @noRd
B2XZG_2d <- function(B, pord = c(2, 2), c = c(10, 10)) {
  c1 <- c[1]
  c2 <- c[2]
  c1c2 <- ncol(B)

  if (c1c2 != c1 * c2) {
    stop("c1 * c2 must me equal to the number of colums of B", call. = FALSE)
  }

  D_1 <- diff(diag(c1), differences = pord[1])
  D_2 <- diff(diag(c2), differences = pord[2])

  P1.svd <- svd(crossprod(D_1))
  P2.svd <- svd(crossprod(D_2))

  U_1s <- (P1.svd$u)[, 1:(c1 - pord[1])] # eigenvectors
  U_1n <- ((P1.svd$u)[, -(1:(c1 - pord[1]))])
  d1 <- (P1.svd$d)[1:(c1 - pord[1])] # eigenvalues

  U_2s <- (P2.svd$u)[, 1:(c2 - pord[2])] # eigenvectors
  U_2n <- ((P2.svd$u)[, -(1:(c2 - pord[2]))])
  d2 <- (P2.svd$d)[1:(c2 - pord[2])] # eigenvalues


  T_n <- kronecker(U_1n, U_2n)

  AUX_1 <- kronecker(U_1n, U_2s)
  AUX_2 <- kronecker(U_1s, U_2n)
  AUX_3 <- kronecker(U_1s, U_2s)

  T_s <- cbind(AUX_1, AUX_2, AUX_3)


  Z <- B %*% T_s
  X <- B %*% T_n

  ####

  d_1s <- diag(P1.svd$d)[1:(c1 - pord[1]), 1:(c1 - pord[1])]
  d_2s <- diag(P2.svd$d)[1:(c2 - pord[2]), 1:(c2 - pord[2])]

  T_1 <- kronecker(diag(pord[1]), d_2s)
  T_2 <- matrix(0, nrow = pord[2] * (c1 - pord[1]), ncol = pord[2] * (c1 - pord[1]))
  T_3 <- kronecker(diag(c1 - pord[1]), d_2s)

  T_21 <- cbind(T_1, matrix(0, nrow = dim(T_1)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_1)[2]))
  T_22 <- cbind(matrix(0, nrow = dim(T_2)[1], ncol = dim(T_1)[2]), T_2, matrix(0, nrow = dim(T_2)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_1)[2] - dim(T_2)[2]))
  T_23 <- cbind(matrix(0, nrow = ((c2 - pord[2]) * (c1 - pord[1])), ncol = (c1 * c2 - pord[1] * pord[2]) - dim(T_3)[2]), T_3)

  H_1 <- matrix(0, nrow = pord[1] * (c2 - pord[2]), ncol = pord[1] * (c2 - pord[2]))
  H_2 <- kronecker(d_1s, diag(pord[2]))
  H_3 <- kronecker(d_1s, diag(c2 - pord[2]))

  H_11 <- cbind(H_1, matrix(0, nrow = dim(H_1)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_1)[2]))
  H_12 <- cbind(matrix(0, nrow = dim(H_2)[1], ncol = dim(H_1)[2]), H_2, matrix(0, nrow = dim(H_2)[1], ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_1)[2] - dim(H_2)[2]))
  H_13 <- cbind(matrix(0, nrow = ((c2 - pord[2]) * (c1 - pord[1])), ncol = (c1 * c2 - pord[1] * pord[2]) - dim(H_3)[2]), H_3)

  L_2 <- rbind(T_21, T_22, T_23)
  L_1 <- rbind(H_11, H_12, H_13)

  t_2 <- diag(L_1)
  t_1 <- diag(L_2)

  G <- list(t_1, t_2)
  names(G) <- c("t_1", "t_2")

  TMatrix <- cbind(T_n, T_s)

  ####

  list(
    X    = X,
    Z    = Z,
    G    = G,
    TMatrix    = TMatrix,
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
