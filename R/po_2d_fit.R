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

    ffpo2d_counter <- 0
    for (ffpo2d_evaluation in evals[grepl("ffpo_2d", names(evals))]) {
      ffpo2d_counter <- ffpo2d_counter + 1

      B_all <- cbind(B_all, ffpo2d_evaluation[["B_ffpo2d"]])
      deglist[[ffpo2d_counter]] <- ffpo2d_evaluation[["nbasis"]][3:4]

      Phi_ffpo[[ffpo2d_counter]] <- ffpo2d_evaluation[["Phi_ffpo2d"]]
      M_ffpo[[ffpo2d_counter]] <- ffpo2d_evaluation[["M_ffpo2d"]] # this is not M, this is M complement
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

  X <- as.matrix(cbind(rep(1, lenght = nrow(X)), X))

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

  # Vp_tmp <- calculate_Vp_tmp(fit, non_special_indices = list(), nf = 0, n_coefs_f = 0)
  #
  # theta_aux <- calculate_theta_aux(fit, non_special_indices = list(), nf = 0, l.f = NULL)
  #
  # covar_theta <- B_res$TMatrix %*% Vp_tmp %*% t(B_res$TMatrix)
  #
  # theta_ffvd <- B_res$TMatrix %*% theta_aux
  #
  # theta_no_functional <- calculate_theta_no_functional(fit, non_special_indices = list())
  #
  # theta_f <- calculate_theta_f(0, l.f = NULL, fit, non_special_indices = list())
  # Beta_ffvd <- calculate_beta_ffvd(nffvd, data, response, M, L_Phi, B_T, theta_ffvd, deglist)
  #
  ffpo_2d_evals <- process_evals(evals, type = "ffpo_2d")

  res <- list(
    fit = fit,
    # Beta = Beta_ffvd,
    theta = theta_ffpo,
    # covar_theta = covar_theta,
    # M = M,
    ffpo_2d_evals = ffpo_2d_evals
    # theta_no_functional = theta_no_functional,
    # theta_f = theta_f
  )

  class(res) <- "po_2d_fit"
  attr(res, "N") <- length(data[[response]])

  res
}

#' @export
summary.po_2d_fit <- function(object, ...) {
  base::summary(object$fit, ...)
}
