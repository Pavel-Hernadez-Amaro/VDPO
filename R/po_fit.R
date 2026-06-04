#' Estimation of functional regression models for partially observed
#' functional data
#'
#' The \code{po_fit} function fits functional regression models for partially
#' observed functional data, where each curve is only observed over part of the
#' common domain.
#'
#' @param formula a formula object with at least one \code{ffpo} term.
#' @param data a \code{list} object containing the response variable
#' and the covariates as the components of the list.
#' @param family a \code{family} object specifying the distribution from which the
#' data originates. The default distribution is \code{\link{gaussian}}.
#' @param offset an offset vector. The default value is \code{NULL}.
#'
#' @return An object of class \code{po_fit}. It is a \code{list} containing the
#' following items:
#'
#' - An item named `fit` of class \code{sop}. See \link[SOP]{sop.fit}.
#' - An item named `Beta` which is a list with one \code{data.frame} per
#' functional term, each containing the grid (`t`), the estimated functional
#' coefficient (`beta`), its standard error (`se`) and the lower and upper
#' limits of the pointwise confidence interval (`lower`, `upper`).
#' - An item named `intercept` which is the estimated intercept of the model.
#' - An item named `theta` which is the basis coefficient vector of the
#' estimated functional coefficient.
#' - An item named `covar_theta` which is the covariance matrix of the basis
#' coefficients, used to build the pointwise confidence intervals.
#' - An item named `M` which holds the observed domain information for each
#' functional term.
#' - An item named `ffpo_evals` which is the result of the evaluations of the
#' `ffpo` terms in the formula.
#'
#' @examples
#' # PARTIALLY OBSERVED FUNCTIONAL DATA EXAMPLE
#'
#' # set seed for reproducibility
#' set.seed(123)
#'
#' # generate example data with partially observed curves
#' sim  <- data_generator_po_1d(n = 50, grid_points = 60)
#' X    <- sim$noisy_curves_miss
#' y    <- sim$response
#' grid <- sim$grid
#' mp   <- sim$missing_points
#'
#' # Fit the model using an 'ffpo' term for the partially observed covariate.
#' # 'nbasis' sets the number of basis functions for the data reconstruction
#' # and for the functional coefficient, respectively.
#' fit  <- po_fit(
#'   response ~ ffpo(X = X, missing_points = mp, grid = grid, nbasis = c(20, 20)),
#'   data = list(response = y, X = X, grid = grid, missing_points = mp)
#' )
#'
#' # Inspect the structure of the returned object
#' str(fit, max.level = 1)
#'
#' # The estimated intercept and the basis coefficients can be accessed directly
#' fit$intercept
#' fit$theta
#'
#' # A summary of the underlying fit can be obtained using the summary function
#' summary(fit)
#'
#' @seealso \code{\link{ffpo}}
#'
#' @export
po_fit <- function(formula, data, family = stats::gaussian(), offset = NULL) {
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

  tf <- stats::terms.formula(formula, specials = c("ffpo"))

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
  nffpo <- sum(grepl("\\bffpo\\(\\b", names(evals)))
  # nffpo_2d <- sum(grepl("\\bffpo_2d\\(\\b", names(evals)))

  if (nffpo == 0) {
    stop(
      "this function should be used with at least one 'ffpo' term",
      call. = FALSE
    )
  }

  X <- c() # matrix
  Z <- c() # matrix
  G <- list() # list

  if (nffpo > 0) {
    B_all <- c()
    deglist <- vector(mode = "list", length = nffpo)
    Phi_ffpo <- vector(mode = "list", length = nffpo) # list of matrices
    M_ffpo <- vector(mode = "list", length = nffpo) # list of matrices
    grid_ffpo <- vector(mode = "list", length = nffpo) # list of grids

    ffpo_counter <- 0
    for (ffpo_evaluation in evals[grepl("ffpo", names(evals))]) {
      ffpo_counter <- ffpo_counter + 1

      B_all <- cbind(B_all, ffpo_evaluation[["B_ffpo"]])
      deglist[[ffpo_counter]] <- ffpo_evaluation[["nbasis"]][2]

      Phi_ffpo[[ffpo_counter]] <- ffpo_evaluation[["Phi"]]
      M_ffpo[[ffpo_counter]] <- ffpo_evaluation[["M"]]
      grid_ffpo[[ffpo_counter]] <- ffpo_evaluation[["grid"]]
    }
    foo <- B2XZG_ffpo(B_all, deglist)

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

  ffpo_evals <- process_evals(evals, type = "ffpo")

  # Reconstruct the functional coefficient(s) and pointwise confidence intervals
  Beta <- calculate_beta_ffpo(theta_ffpo, covar_theta, Phi_ffpo, grid_ffpo, deglist)

  res <- list(
    fit = fit,
    Beta = Beta,
    intercept = intercept,
    theta = theta_ffpo,
    covar_theta = covar_theta,
    M = M_ffpo,
    ffpo_evals = ffpo_evals
  )

  class(res) <- "po_fit"
  attr(res, "N") <- length(data[[response]])

  res
}

#' Reconstruct the functional coefficient and pointwise confidence intervals for
#' \code{ffpo} terms
#'
#' @noRd
calculate_beta_ffpo <- function(theta, covar_theta, Phi_list, grid_list, deglist,
                                level = 0.95) {
  nterms <- length(Phi_list)
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- vector("list", nterms)
  start <- 1
  for (j in seq_len(nterms)) {
    cj <- deglist[[j]]
    idx <- start:(start + cj - 1)
    Phi_B <- Phi_list[[j]]$B
    theta_j <- theta[idx]
    covar_j <- covar_theta[idx, idx, drop = FALSE]
    beta <- as.numeric(Phi_B %*% theta_j)
    var_beta <- rowSums((Phi_B %*% covar_j) * Phi_B)
    se <- sqrt(pmax(var_beta, 0))
    out[[j]] <- data.frame(
      t = grid_list[[j]],
      beta = beta,
      se = se,
      lower = beta - z * se,
      upper = beta + z * se
    )
    start <- start + cj
  }
  out
}

#' @export
summary.po_fit <- function(object, ...) {
  base::summary(object$fit, ...)
}
