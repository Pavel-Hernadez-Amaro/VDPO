sum((data$Beta[,,1] - res$Beta_ffvd[[1]])^2, na.rm = TRUE)/(100*101)
sum((data$y - res$fit$fitted.values)^2, na.rm = TRUE)/(100)

plotly::plot_ly(z = data$Beta[,,2], type = "surface")
plotly::plot_ly(z = res$Beta_ffvd[[1]], type = "surface")




ffvd_old <- function(X, grid, nbasis = c(30, 50, 30), bdeg = c(3, 3, 3)) {

  if (missing(grid)) {
    grid <- seq_len(ncol(X))
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

  X_hat=matrix(nrow = N, ncol = ncol(X))

  c1 <- nbasis[1]
  c2 <- nbasis[2]
  c3 <- nbasis[3]

  K <- NULL
  rng <- cbind(grid[M[,1]],grid[M[,2]])

  L_Phi <- vector(mode = "list", length = N)
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

    X_hat[i, M[i,1]:M[i,2]]=aux_3$fit$fitted.values

    A[i, ((c1 * (i - 1)) + 1):(i * c1)] <- aux_3$theta

    ############### HERE WE CREATE THE MARGINAL BASIS FOR THE t VARIABLE In B(t,T)

    c_t <- c2 - bdeg[2] # EQUAL TO THE NUMBER OF INNER KNOTS + 1

    L_Phi[[i]] <- bspline(grid[M[i, 1]:M[i, 2]], XL, XR, c_t, bdeg[2])
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

  # L_X_all <- bspline(grid, min(grid) - 1e-6, max(grid) + 1e-6, c_t, bdeg[1])

  # PERFORMING THE INNER PRODUCT


  # need to rewrite this for statement
  for (i in 1:N) {
    PROD <- partial_inprod(
      n_intervals   = sub,
      knots1        = L_X[[i]]$knots,
      knots2        = L_Phi[[i]]$knots,
      bdeg          = bdeg[1:2],
      spline_domain = B_T$B[i, , drop = FALSE],
      rng           = c(grid[M[i, 1]], grid[M[i, 2]])
    )
    PROD <- PROD / grid[(M[i, 2] - M[i, 1] + 1)]

    K <- rbind(K, PROD)
  }

  B <- A %*% K

  list(
    B = B,
    X_hat  = X_hat,
    L_Phi  = L_Phi,
    B_T    = B_T,
    M      = M,
    nbasis = nbasis
  )
}


VDPO <- function(formula, data, family = stats::gaussian(), offset = NULL) {
  if (inherits(formula, "character")) {
    formula <- stats::as.formula(formula)
  }
  if (inherits(data, what = "data.frame")) {
    data <- as.data.frame(data)
  } else {
    stop("The data specified in the 'data' argument should be a data frame",
         call. = FALSE
    )
  }

  na_indices <- c()
  for (name in names(data)) {
    if (is.null(dim(data[[name]]))) {
      na_indices <- c(na_indices, which(is.na(data[[name]])))
    }
  }

  na_indices <- unique(na_indices)

  if (length(na_indices) != 0) {
    data <- data[-na_indices, ]
  }

  vdpoenv <- environment(formula)
  vdpons <- loadNamespace("VDPO")

  for (var in names(data)) {
    vdpoenv[[var]] <- data[[var]]
  }
  nobs <- nrow(data)

  if (is.null(offset)) {
    offset <- rep(0L, nobs)
  }

  tf <- stats::terms.formula(formula, specials = c("ffvd_old", "ffpo", "ffpo_2d", "f"))

  terms <- attr(tf, "term.labels")
  nterms <- length(terms)
  specials_indices <- attr(tf, "specials") # indices for the special terms

  if (attr(tf, "response")) {
    response <- as.character(attr(tf, "variables"))[2]
    variables <- as.character(attr(tf, "variables"))[-c(1, 2)]

    # If the response exists, we need to lower in the index for every
    # special term

    specials_indices <- lapply(
      specials_indices,
      function(x) if (is.null(x)) NA else x - 1
    )
  } else {
    variables <- as.character(attr(tf, "variables"))[-1]

    # No need to lower the index of the specials terms if the response exists
    specials_indices <- lapply(
      specials_indices,
      function(x) if (is.null(x)) NA else x
    )
  }
  all_indices <- seq_along(variables)
  non_special_indices <- setdiff(all_indices, unlist(specials_indices))

  specials <- length(unlist(stats::na.omit(specials_indices))) > 0

  # if (specials) {
  #   non_funcional_variables <- variables[-specials_indices]
  #   functional_variables    <- variables[specials_indices]
  # } else {
  #   # pass
  # }

  evals <- lapply(
    terms,
    function(term) eval(parse(text = term), envir = vdpoenv, enclos = vdpons)
  )
  names(evals) <- terms

  nf       <- sum(grepl("\\bf\\(\\b", names(evals)))
  nffvd    <- sum(grepl("\\bffvd_old\\(\\b", names(evals)))
  nffpo    <- sum(grepl("\\bffpo\\(\\b", names(evals)))
  nffpo_2d <- sum(grepl("\\bffpo_2d\\(\\b", names(evals)))

  ## TODO add the f function to the package namespace

  if (nffvd == 0 && nffpo == 0 && nffpo_2d == 0) {
    stop("this function should be used with at least one 'ffvd_old', 'ffpo' or 'ffpo_2d' term",
         call. = FALSE
    )
  }

  X <- c() # matrix
  Z <- c() # matrix
  G <- list() # list

  l.f <- NULL
  if (nf > 0) {
    names.fun <- NULL
    for (f_index in which(grepl("\\bf\\(\\b", names(evals)))) {
      l.f[[f_index]] <- evals[[f_index]]
      names.fun <- c(names.fun, terms[[f_index]])
      form <- stats::as.formula(paste("~", terms[[f_index]], sep = ""))
      aa <- switch(as.character(l.f[[f_index]]$dim),
                   `3` = {
                     l.f[[f_index]]$Xmat <- construct.3D.pspline(
                       form,
                       data
                     )
                   },
                   `2` = {
                     l.f[[f_index]]$Xmat <- construct.2D.pspline(
                       form,
                       data
                     )
                   },
                   `1` = {
                     l.f[[f_index]]$Xmat <- construct.1D.pspline(
                       form,
                       data
                     )
                   }
      )
      X <- cbind(X, l.f[[f_index]]$Xmat$X)
      Z <- cbind(Z, l.f[[f_index]]$Xmat$Z)
      G <- c(G, list(l.f[[f_index]]$Xmat$g))
    }
    l.f <- l.f[lengths(l.f) != 0]

    if (length(G) == 1 && length(G[[1]]) == 1) {
      # G <- list(unlist(G))
      G <- unname(G[[1]])
      # names(G) <- names(G[[1]])
    } else {
      G <- construct.capital.lambda(G)
    }

    aux_sum <- 0
    for (i in seq_along(l.f)) {
      aux_prod <- 1
      for (j in seq_along(l.f[[i]]$nseg)) {
        aux_prod <- aux_prod * (l.f[[i]]$nseg[j] + l.f[[i]]$degree[j])
      }
      aux_sum <- aux_sum + aux_prod
    }
  }

  if (nffvd > 0) {
    B_all   <- c()
    deglist <- vector(mode = "list", length = nffvd)
    L_Phi   <- vector(mode = "list", length = nffvd) # list of matrices
    B_T     <- vector(mode = "list", length = nffvd) # list of matrices
    M       <- vector(mode = "list", length = nffvd) # list of matrices
    # TMatrix <- vector(mode = "list", length = nffvd) # matrix

    ffvd_counter <- 0
    for (ffvd_evaluation in evals[grepl("ffvd_old", names(evals))]) {
      ffvd_counter <- ffvd_counter + 1

      B_all <- cbind(B_all, ffvd_evaluation[["B"]])
      deglist[[ffvd_counter]] <- ffvd_evaluation[["nbasis"]][2:3]

      L_Phi[[ffvd_counter]] <- ffvd_evaluation[["L_Phi"]]
      B_T[[ffvd_counter]] <- ffvd_evaluation[["B_T"]]
      M[[ffvd_counter]] <- ffvd_evaluation[["M"]]
      # TMatrix[[ffvd_counter]] <- ffvd_evaluation[["TMatrix"]]
    }
    foo <- B2XZG(B_all, deglist)

    X <- cbind(X, foo$X_ffvd)
    Z <- cbind(Z, foo$Z_ffvd)
    G <- c(G, foo$G_ffvd)
  }

  if (nffpo > 0) {
    B_all <- c()
    deglist <- vector(mode = "list", length = nffpo)
    Phi_ffpo <- vector(mode = "list", length = nffpo) # list of matrices
    M_ffpo <- vector(mode = "list", length = nffpo) # list of matrices
    # TMatrix <- vector(mode = "list", length = nffvd) # matrix

    ffpo_counter <- 0
    for (ffpo_evaluation in evals[grepl("ffpo", names(evals))]) {
      ffpo_counter <- ffpo_counter + 1

      B_all <- cbind(B_all, ffpo_evaluation[["B_ffpo"]])
      deglist[[ffpo_counter]] <- ffpo_evaluation[["nbasis"]][1]

      Phi_ffpo[[ffpo_counter]] <- ffpo_evaluation[["Phi"]]
      M_ffpo[[ffpo_counter]] <- ffpo_evaluation[["M"]]
    }
    foo <- B2XZG_ffpo(B_all, deglist)

    X <- cbind(X, foo$X)
    Z <- cbind(Z, foo$Z)
    G <- c(G, foo$G)
  }

  if (nffpo_2d > 0) {
    B_all <- c()
    deglist <- vector(mode = "list", length = nffpo_2d)
    Phi_ffpo <- vector(mode = "list", length = nffpo_2d) # list of matrices
    M_ffpo <- vector(mode = "list", length = nffpo_2d) # list of matrices
    # TMatrix <- vector(mode = "list", length = nffvd) # matrix

    ffpo2d_counter <- 0
    for (ffpo2d_evaluation in evals[grepl("ffpo_2d", names(evals))]) {
      ffpo2d_counter <- ffpo2d_counter + 1

      B_all <- cbind(B_all, ffpo2d_evaluation[["B_ffpo2d"]])
      deglist[[ffpo2d_counter]] <- ffpo2d_evaluation[["nbasis"]][3:4]

      Phi_ffpo[[ffpo2d_counter]] <- ffpo2d_evaluation[["Phi_ffpo2d"]]
      M_ffpo[[ffpo2d_counter]] <- ffpo2d_evaluation[["M_ffpo2d"]] # this is not M, this is M complement
    }
    foo <- B2XZG_2d(B_all, pord = c(2, 2), c = deglist[[1]]) # we need to fix this to work like deglist in ffvd_old || multivariate case

    X <- cbind(X, foo$X)
    Z <- cbind(Z, foo$Z)
    G <- c(G, foo$G)
  }

  for (column_index in rev(non_special_indices)) {
    X <- cbind(evals[[column_index]], X)
  }

  if (is.null(Z)) {
    stop("check the 'lm' function from the 'stats' package", call. = TRUE)
  }

  X <- as.matrix(cbind(rep(1, lenght = nrow(X)), X))


  # We need to add zeros to the right of the f Gs and to the left of the ffvd_old

  if (nf > 0) {
    for (i in 1:nf) {
      G[[i]] <- c(G[[i]], rep(0, ncol(Z) - length(G[[i]])))
    }
  }
  # Adding zeros to the left for the ffvd_old
  # for (i in (nf+1):(2*(nf+nffvd)-1)) {
  #   G[[i]] <- c(rep(0, ncol(Z) - length(G[[i]])), G[[i]])
  # }

  fit <- sop.fit(
    X = X, Z = Z, G = G,
    y = data[[response]], family = family,
    control = list(trace = FALSE), offset = offset
  )

  intercept <- fit$b.fixed[1]
  b_fixed_tmp <- fit$b.fixed[-1]
  Vp_tmp <- fit$Vp[-1, -1]

  if (length(non_special_indices) == 0 && nf == 0 && nffvd > 0 && nffpo == 0 && nffpo_2d == 0) {
    # only ffvd_old

    theta_aux <- c(b_fixed_tmp, fit$b.random)
    covar_theta <- foo$TMatrix %*% Vp_tmp %*% t(foo$TMatrix)
    std_error_theta <- sqrt(diag(foo$TMatrix %*% Vp_tmp %*% t(foo$TMatrix)))
  } else if (length(non_special_indices) == 0 && nf == 0 && nffvd == 0 && nffpo > 0 && nffpo_2d == 0) {
    # only ffpo

    theta_aux <- c(b_fixed_tmp, fit$b.random)
  } else if (length(non_special_indices) == 0 && nf == 0 && nffvd == 0 && nffpo == 0 && nffpo_2d > 0) {
    # only ffpo_2d

    theta_aux <- c(b_fixed_tmp, fit$b.random)
  }
  # return(list(
  #   TMatrix = foo$TMatrix, theta_aux = theta_aux
  # ))
  theta <- foo$TMatrix %*% theta_aux
  if (nffvd > 0) {
    Beta_ffvd <- lapply(
      M,
      function(x) {
        matrix(NA,
               nrow = length(data[[response]]),
               ncol = max(x)
        )
      }
    )

    it <- 1
    for (j in 1:nffvd) {
      for (i in 1:length(data[[response]])) {
        Beta_ffvd[[j]][i, M[[j]][i, 1]:M[[j]][i, 2]] <-
          as.matrix(kronecker(L_Phi[[j]][[i]]$B[M[[j]][i, 1]:M[[j]][i, 2],], t(B_T[[j]]$B[i, ]))) %*% theta[it:(it + prod(deglist[[j]]) - 1)]
      }
      it <- it + prod(deglist[[j]])
    }
  }

  if (length(non_special_indices) == 0 && nf == 0 && nffvd > 0 && nffpo == 0 && nffpo_2d == 0) {
    # only ffvd_old

    res <- list(
      fit              = fit,
      theta_ffvd       = theta,
      covar_theta      = covar_theta,
      std_error_theta  = std_error_theta,
      Beta_ffvd        = Beta_ffvd,
      M_ffvd           = M
    )

  } else if (length(non_special_indices) == 0 && nf == 0 && nffvd == 0 && nffpo > 0 && nffpo_2d == 0) {
    # only ffpo

    res <- list(
      fit              = fit,
      theta_ffpo       = theta
      # Beta_ffpo        = Beta_ffpo,
      # M_ffpo           = M_ffpo
    )

  } else if (length(non_special_indices) == 0 && nf == 0 && nffvd == 0 && nffpo == 0 && nffpo_2d > 0) {
    # only ffpo_2d

    res <- list(
      fit              = fit,
      theta_ffpo2d     = theta
      # Beta_ffpo2d      = Beta_ffvd,
      # M_ffpo2d         = M
    )
  }


  class(res) <- "VDPO"
  attr(res, "N") <- length(data[[response]])
  res
}

#' @export
plot.VDPO <- function(x, beta_index = 1, ...) {
  Beta <- NULL

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2' is required for this functionality", call. = FALSE)
  }

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("package 'RColorBrewer' is required for this functionality", call. = FALSE)
  }

  if (beta_index < 1 || beta_index > length(x$M_ffvd)) {
    stop(
      "'beta_index' should be between 1 and the number of variable domain
         functional variables used in the formula",
      call. = FALSE
    )
  }

  N <- attr(x, "N")

  max_M <- max(x$M_ffvd[[beta_index]])
  T_dat <- IND <- NULL
  t_dat <- rep(1:max_M, N)

  Beta_estimated <- t(x$Beta_ffvd[[beta_index]])
  dim(Beta_estimated) <- c(nrow(x$Beta_ffvd[[beta_index]]) * ncol(x$Beta_ffvd[[beta_index]]), 1)

  for (ind in 1:N) {
    T_dat <- c(T_dat, rep(x$M_ffvd[[beta_index]][ind, 2], max_M))
    IND <- c(IND, rep(ind, x$M_ffvd[[beta_index]][ind, 2]))
  }

  Heat_map_data <- data.frame(t = t_dat, M = T_dat, Beta = Beta_estimated)
  Heat_map_data <- Heat_map_data[t_dat <= T_dat, ]
  Heat_map_data[["IND"]] <- IND

  lims <- range(Heat_map_data$Beta)

  ggplot2::ggplot(Heat_map_data, ggplot2::aes(x = t, y = IND)) +
    ggplot2::geom_tile(ggplot2::aes(colour = Beta, fill = Beta)) +
    ggplot2::scale_fill_gradientn(
      name     = "",
      limits   = lims,
      colours  = rev(RColorBrewer::brewer.pal(11, "Spectral")),
      na.value = "white"
    ) +
    ggplot2::scale_colour_gradientn(
      name    = "",
      limits  = lims,
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "T") +
    ggplot2::ggtitle(paste0("FFVD HeatMap (Beta ", beta_index, ")"))
}
