#' Estimation of functional regression models for partially observed
#' bidimensional functional data
#'
#' The \code{po_2d_fit} function fits functional regression models for partially
#' observed bidimensional functional data, where each surface is only observed
#' over part of the common domain.
#'
#' @param formula a formula object with at least one \code{ffpo_2d} term.
#' @param data a \code{list} object containing the response variable
#' and the covariates as the components of the list.
#' @param family a \code{family} object specifying the distribution from which the
#' data originates. The default distribution is \code{\link{gaussian}}.
#' @param offset an offset vector. The default value is \code{NULL}.
#'
#' @return An object of class \code{po_2d_fit}. It is a \code{list} containing the
#' following items:
#'
#' - An item named `fit` of class \code{sop}. See \link[SOP]{sop.fit}.
#' - An item named `theta` which is the basis coefficient vector of the
#' estimated bidimensional functional coefficient.
#' - An item named `ffpo_2d_evals` which is the result of the evaluations of the
#' `ffpo_2d` terms in the formula.
#'
#' @examples
#' # PARTIALLY OBSERVED BIDIMENSIONAL FUNCTIONAL DATA EXAMPLE
#' \donttest{
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # generate example data with partially observed surfaces
#' sim  <- data_generator_po_2d(n = 30, grid_x = 8, grid_y = 8)
#' X    <- sim$noisy_surfaces_miss
#' y    <- sim$response
#' mp   <- sim$missing_points
#' mpts <- sim$miss_points
#'
#' # Fit the model using an 'ffpo_2d' term for the partially observed surfaces.
#' # 'miss_points' and 'missing_points' describe the missing observations in the
#' # two complementary formats expected by 'ffpo_2d'.
#' fit  <- po_2d_fit(
#'   response ~ ffpo_2d(
#'     X = X, miss_points = mpts, missing_points = mp, nbasis = rep(5, 4)
#'   ),
#'   data = list(response = y, X = X)
#' )
#'
#' # Inspect the structure of the returned object
#' str(fit, max.level = 1)
#'
#' # The basis coefficients of the functional coefficient can be accessed directly
#' fit$theta
#'
#' # A summary of the underlying fit can be obtained using the summary function
#' summary(fit)
#' }
#'
#' @seealso \code{\link{ffpo_2d}}
#'
#' @export
po_2d_fit <- function(formula, data, family = stats::gaussian(), offset = NULL) {
  if (inherits(formula, "character")) {
    formula <- stats::as.formula(formula)
  }

  if (!inherits(data, what = "list")) {
    stop("the 'data' argument should be a list", call. = FALSE)
  }

  if (length(data) < 2) {
    stop("'data' should have at least a response and a covariate", call. = FALSE)
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
  vdpoenv$f <- SOP::f
  vdpons <- loadNamespace("VDPO")

  for (var in names(data)) {
    vdpoenv[[var]] <- data[[var]]
  }

  tf <- stats::terms.formula(formula, specials = c("ffpo_2d"))

  terms <- attr(tf, "term.labels")
  nterms <- length(terms)
  specials_indices <- attr(tf, "specials") # indices for the special terms

  if (attr(tf, "response")) {
    response <- as.character(attr(tf, "variables"))[2]
    variables <- as.character(attr(tf, "variables"))[-c(1, 2)]

    specials_indices <- lapply(
      specials_indices,
      function(x) if (is.null(x)) NA else x - 1
    )
  } else {
    variables <- as.character(attr(tf, "variables"))[-1]
    specials_indices <- lapply(
      specials_indices,
      function(x) if (is.null(x)) NA else x
    )
  }

  all_indices <- seq_along(variables)
  non_special_indices <- setdiff(all_indices, unlist(specials_indices))

  specials <- length(unlist(stats::na.omit(specials_indices))) > 0

  evals <- lapply(
    terms,
    function(term) eval(parse(text = term), envir = vdpoenv, enclos = vdpons)
  )
  names(evals) <- terms

  # nf <- sum(grepl("\\bf\\(\\b", names(evals)))
  # nffpo    <- sum(grepl("\\bffpo\\(\\b", names(evals)))
  nffpo_2d <- sum(grepl("\\bffpo_2d\\(\\b", names(evals)))

  if (nffpo_2d == 0) {
    stop(
      "this function should be used with at least one 'ffpo_2d' term",
      call. = FALSE
    )
  }

  X <- c() # matrix
  Z <- c() # matrix
  G <- list() # list

  if (nffpo_2d > 0) {
    B_all <- c()
    deglist <- vector(mode = "list", length = nffpo_2d)
    Phi_ffpo <- vector(mode = "list", length = nffpo_2d) # list of matrices
    M_ffpo <- vector(mode = "list", length = nffpo_2d) # list of matrices
    grid_ffpo <- vector(mode = "list", length = nffpo_2d) # list of x/y grids

    ffpo2d_counter <- 0
    for (ffpo2d_evaluation in evals[grepl("ffpo_2d", names(evals))]) {
      ffpo2d_counter <- ffpo2d_counter + 1

      B_all <- cbind(B_all, ffpo2d_evaluation[["B_ffpo2d"]])
      deglist[[ffpo2d_counter]] <- ffpo2d_evaluation[["nbasis"]][3:4]

      Phi_ffpo[[ffpo2d_counter]] <- ffpo2d_evaluation[["Phi_ffpo2d"]]
      M_ffpo[[ffpo2d_counter]] <- ffpo2d_evaluation[["M_ffpo2d"]] # this is not M, this is M complement
      grid_ffpo[[ffpo2d_counter]] <- list(
        x = ffpo2d_evaluation[["points_x"]],
        y = ffpo2d_evaluation[["points_y"]]
      )
    }
    foo <- B2XZG_2d(B_all, pord = c(2, 2), c = deglist[[1]]) # we need to fix this to work like deglist in ffvd || multivariate case

    X <- cbind(X, foo$X)
    Z <- cbind(Z, foo$Z)
    G <- c(G, foo$G)
  }

  for (column_index in rev(non_special_indices)) {
    X <- cbind(evals[[column_index]], X)
  }

  if (is.null(Z)) {
    stop("please use the 'lm' function from the 'stats' package", call. = TRUE)
  }

  X <- as.matrix(cbind(rep(1, length.out = nrow(X)), X))

  if (is.null(offset)) {
    nobs <- length(data[[response]])
    offset <- rep(0L, nobs)
  }

  fit <- sop.fit(
    X = X,
    Z = Z,
    G = G,
    y = data[[response]],
    family = family,
    control = list(trace = FALSE),
    offset = offset
  )

  b_fixed_tmp <- fit$b.fixed[-1]
  theta_aux <- c(b_fixed_tmp, fit$b.random)
  theta_ffpo <- foo$TMatrix %*% theta_aux

  intercept <- fit$b.fixed[1]

  # Covariance of the spline coefficients (drop the intercept column/row, as in vd_fit)
  Vp_tmp <- fit$Vp[-1, -1]
  covar_theta <- foo$TMatrix %*% Vp_tmp %*% t(foo$TMatrix)

  ffpo_2d_evals <- process_evals(evals, type = "ffpo_2d")

  # Reconstruct the coefficient surface(s) and pointwise confidence intervals
  Beta <- calculate_beta_ffpo_2d(theta_ffpo, covar_theta, Phi_ffpo, grid_ffpo, deglist)

  res <- list(
    fit = fit,
    Beta = Beta,
    intercept = intercept,
    theta = theta_ffpo,
    covar_theta = covar_theta,
    M = M_ffpo,
    ffpo_2d_evals = ffpo_2d_evals
  )

  class(res) <- "po_2d_fit"
  attr(res, "N") <- length(data[[response]])

  res
}

#' @export
summary.po_2d_fit <- function(object, ...) {
  base::summary(object$fit, ...)
}

#' Reconstruct the coefficient surface and pointwise confidence intervals for
#' \code{ffpo_2d} terms
#'
#' @noRd
calculate_beta_ffpo_2d <- function(theta, covar_theta, Phi_list, grid_list, deglist,
                                   level = 0.95) {
  nterms <- length(Phi_list)
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- vector("list", nterms)
  start <- 1
  for (j in seq_len(nterms)) {
    cj <- prod(deglist[[j]])
    idx <- start:(start + cj - 1)
    Phi_B <- Phi_list[[j]]
    theta_j <- theta[idx]
    covar_j <- covar_theta[idx, idx, drop = FALSE]
    beta_vec <- as.numeric(Phi_B %*% theta_j)
    var_vec <- rowSums((Phi_B %*% covar_j) * Phi_B)
    se_vec <- sqrt(pmax(var_vec, 0))
    nx <- length(grid_list[[j]]$x)
    ny <- length(grid_list[[j]]$y)
    out[[j]] <- list(
      x = grid_list[[j]]$x,
      y = grid_list[[j]]$y,
      beta = matrix(beta_vec, nrow = nx, ncol = ny),
      se = matrix(se_vec, nrow = nx, ncol = ny),
      lower = matrix(beta_vec - z * se_vec, nrow = nx, ncol = ny),
      upper = matrix(beta_vec + z * se_vec, nrow = nx, ncol = ny)
    )
    start <- start + cj
  }
  out
}
