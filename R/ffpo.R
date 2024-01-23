#' ffpo
#'
#' @param X .
#' @param grid .
#' @param bidimensional_grid .
#' @param nbasis .
#' @param bdeg .
#'
#' @return .
#'
#' @export
ffpo <- function(X, grid, bidimensional_grid = FALSE, nbasis = c(30, 30), bdeg = c(3, 3)) {
  if (!is.matrix(X)) {
    stop("argument 'X' should be a matrix", call. = FALSE)
  }

  if (!bidimensional_grid) {
    grid <- c(grid)
    grid <- grid[!is.na(grid)]
  }

  if (!is.null(dim(grid)) && length(dim(grid)) > 2) {
    stop("the 'grid' parameter should be at most 2-dimensional", call. = FALSE)
  }

  grid <- as.matrix(grid)
  is_grid_matrix <- length(setdiff(dim(grid),1)) == 2

  SUB <- 500
  K   <- NULL
  N   <- nrow(X)
  c1  <- nbasis[1]
  c2  <- nbasis[2]
  A <- matrix(0, nrow = N, ncol = N*c1)

  if (!is_grid_matrix) {
    M <- t(apply(X, 1, function(x) range(which(!is.na(x)))))
  } else {
    grid_all <- sort(unique(c(t(grid)))) # Mathematical union of the grid rows
    M <- c(apply(grid, 1, function(x) sum(!is.na(x))))
  }

  rng <- matrix(0, ncol = 2, nrow = N)

  L_X       <- vector(mode = "list", length = N)
  L_y       <- vector(mode = "list", length = N)
  L_theta   <- vector(mode = "list", length = N)

  for (i in 1:N) {
    if (!is_grid_matrix) {
      rng[i,] <- c(grid[M[i,1]], grid[M[i,2]]) #c(range(grid)[1],range(grid)[2])
    } else {
      rng[i,] <- c(grid[i,1], grid[i,M[i]])
    }

    XL <- rng[i,1] - 1e-06
    XR <- rng[i,2] + 1e-06

    c <- c1 - bdeg[1] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    if (!is_grid_matrix) {
      L_X[[i]] <- bspline(grid[M[i,1]:M[i,2]], XL, XR, c, bdeg[1])
      response_x <- X[i,M[i,1]:M[i,2]]
    } else {
      L_X[[i]] <- bspline(grid[i, 1:M[i]], XL, XR, c, bdeg[1])
      response_x <- X[i, 1:M[i]]
    }

    # Estimating the data coefficients (Matrix A)

    aux <- L_X[[i]]$B
    aux_2 <- B2XZG_1d(aux,2,c1)
    aux_3 <- XZG2theta_1d(X = aux_2$X, Z = aux_2$Z, G = aux_2$G, TMatrix = aux_2$T, y = response_x)

    A[i,((c1*(i-1))+1):(i*c1)] <- aux_3$theta

    L_theta[[i]] <- aux_3$theta
    L_y[[i]] <- L_X[[i]]$B%*%L_theta[[i]]
  }


  # The knots for the functional coefficient must be maximum
  rng_t <- if (!is_grid_matrix) range(grid) else range(grid_all, na.rm = TRUE)
  XL_t <- rng_t[1] - 1e-06
  XR_t <- rng_t[2] + 1e-06
  c_t <- c2 - bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1
  Phi <- bspline(if (!is_grid_matrix) grid else grid_all,
                 XL_t, XR_t, c_t, bdeg[2])

  res <- matrix(nrow = N, ncol = c2)

  for (i in 1:N) {
    aux_knots <- if (!is_grid_matrix) M[i,1]:M[i,2] else which(grid_all %in% grid[i, ])
    Phi_short <- Phi$B[aux_knots,]

    aux_del <- which(colSums(Phi_short) == 0)

    if (length(aux_del) != 0) {
      Phi_short <- as.matrix(Phi_short[,-aux_del])
    }
    aux_range <- if (!is_grid_matrix) M[i, 1]:M[i, 2] else i
    min_knot <- max(which(Phi$knots <= range(grid[aux_range, ],na.rm = TRUE)[1])) - 3
    max_knot <- min(which(Phi$knots >= range(grid[aux_range, ],na.rm = TRUE)[2])) + 3

    if (min_knot<1) {
      min_knot <- 1
    }

    if (max_knot > length(Phi$knots)) {
      max_knot <- length(Phi$knots)
    }

    new_knots <- Phi$knots[min_knot:max_knot]

    B_X_a   <- splines::spline.des(L_X[[i]]$knots, rng[i,1], bdeg[1]+1, 0*rng[i,1])$design
    B_Phi_a <- splines::spline.des(new_knots, rng[i,1], bdeg[2]+1, 0*rng[i,1])$design

    B_X_b   <- splines::spline.des(L_X[[i]]$knots, rng[i,2], bdeg[1]+1, 0*rng[i,2])$design
    B_Phi_b <- splines::spline.des(new_knots, rng[i,2], bdeg[2]+1, 0*rng[i,2])$design

    n <- 2*SUB # Intervals for the Simpson method

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

      if (j %% 2 != 0) {
        XI2 <- XI2 + Fx
      } else {
        XI1 <- XI1 + Fx
      }
    }

    XI <- matrix(0, nrow = c1, ncol = c2)

    if (length(aux_del) != 0) {
      XI[,-aux_del] <- width * (XIa + XIb + 2*XI2 + 4*XI1)/3
    } else {
      XI <- width * (XIa + XIb + 2*XI2 + 4*XI1)/3
    }

    res[i,] <- t(L_theta[[i]]) %*% XI
  }

  list(
    B_ffpo = res,
    A_ffpo = A,
    K_ffpo = K,
    nbasis = nbasis,
    Phi    = Phi,
    M      = M
  )
}
