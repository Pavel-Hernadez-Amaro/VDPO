#' Estimation of the generalized additive functional regression models for
#' variable domain functional data
#'
#' The \code{vd_fit} function fits generalized additive functional regression models
#' for variable domain functional data.
#'
#' @param formula a formula object with at least one \code{ffvd} term.
#' @param data a \code{list} object containing the response variable
#' and the covariates as the components of the list.
#' @param family a \code{family} object specifying the distribution from which the
#' data originates. The default distribution is \code{\link{gaussian}}.
#' @param offset An offset vector. The default value is \code{NULL}.
#'
#' @return An object of class \code{vd_fit}. It is a \code{list} containing the following items:
#'
#' - An item named `fit` of class \code{sop}. See \link[SOP]{sop.fit}.
#' - An item named `Beta` which is the estimated functional coefficient.
#' - An item named `theta` which is the basis coefficient of `Beta`.
#' - An item named `covar_theta` which is the covariance matrix of `theta`.
#' - An item named `M` which is the number of observations points for each curve.
#' - An item named `ffvd_evals` which is the result of the evaluations of the `ffvd`
#' terms in the formula.
#
#' @examples
#' # VARIABLE DOMAIN FUNCTIONAL DATA EXAMPLE
#'
#' # set seed for reproducibility
#' set.seed(42)
#'
#' # generate example data
#' data <- data_generator_vd(
#'   N = 100,
#'   J = 100,
#'   beta_index = 1,
#'   use_x = TRUE,
#'   use_f = TRUE,
#' )
#'
#' # Define a formula object that specifies the model behavior.
#' # The formula includes a functional form of the variable 'X_se' using 'ffvd'
#' # with a non-default number of basis functions ('nbasis' is set to c(10, 10, 10)).
#' # Additionally, it includes a smooth function 'f' applied to 'x2' with 10 segments ('nseg = 10'),
#' # a second-order penalty ('pord = 2'), and cubic splines ('degree = 3').
#' # The model also contains the linear term 'x1'.
#' formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10)) + f(x2, nseg = 10, pord = 2, degree = 3) + x1
#'
#' # We can fit the model using the data and the formula
#' res <- vd_fit(formula = formula, data = data)
#'
#' # Some important parameters of the model can be accesed as follows
#' res$Beta # variable domain functional coefficient
#' res$fit$fitted.values # estimated response variable
#'
#' # Also, a summary of the fit can be accesed using the summary function
#' summary(res)
#'
#' # And a heatmap for an specific beta can be obtained using the plot function
#' plot(res, beta_index = 1)
#'
#' @seealso \code{\link{ffvd}}
#'
#' @export
vd_fit <- function(formula, data, family = stats::gaussian(), offset = NULL) {
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

  tf <- stats::terms.formula(formula, specials = c("ffvd", "f"))

  terms <- attr(tf, "term.labels")
  nterms <- length(terms)
  specials_indices <- attr(tf, "specials") # indices for the special terms

  if (attr(tf, "response")) {
    # > formula <- y ~ x + z + 1
    # > "list" "y"    "x"    "z"
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

  evals <- lapply(
    terms,
    function(term) eval(parse(text = term), envir = vdpoenv, enclos = vdpons)
  )
  names(evals) <- terms

  nf <- sum(grepl("\\bf\\(\\b", names(evals)))
  nffvd <- sum(grepl("\\bffvd\\(\\b", names(evals)))

  if (nffvd == 0) {
    stop("this function should be used with at least one 'ffvd' term",
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
            list_to_df(data, response)
          )
        },
        `2` = {
          l.f[[f_index]]$Xmat <- construct.2D.pspline(
            form,
            list_to_df(data, response)
          )
        },
        `1` = {
          l.f[[f_index]]$Xmat <- construct.1D.pspline(
            form,
            list_to_df(data, response)
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
      G <- unname(construct.capital.lambda(G))
    }

    n_coefs_f <- 0
    for (i in seq_along(l.f)) {
      aux_prod <- 1
      for (j in seq_along(l.f[[i]]$nseg)) {
        aux_prod <- aux_prod * (l.f[[i]]$nseg[j] + l.f[[i]]$degree[j] - 1)
      }
      n_coefs_f <- n_coefs_f + aux_prod
    }
  }

  if (nffvd > 0) {
    B_all <- c()
    deglist <- vector(mode = "list", length = nffvd)
    L_Phi <- vector(mode = "list", length = nffvd) # list of matrices
    B_T <- vector(mode = "list", length = nffvd) # list of matrices
    M <- vector(mode = "list", length = nffvd) # list of matrices

    ffvd_counter <- 0
    for (ffvd_evaluation in evals[grepl("ffvd", names(evals))]) {
      ffvd_counter <- ffvd_counter + 1

      B_all <- cbind(B_all, ffvd_evaluation[["B"]])
      deglist[[ffvd_counter]] <- ffvd_evaluation[["nbasis"]][2:3]

      L_Phi[[ffvd_counter]] <- ffvd_evaluation[["L_Phi"]]
      B_T[[ffvd_counter]] <- ffvd_evaluation[["B_T"]]
      M[[ffvd_counter]] <- ffvd_evaluation[["M"]]
      # TMatrix[[ffvd_counter]] <- ffvd_evaluation[["TMatrix"]]
    }

    B_res <- B2XZG(B_all, deglist)

    X <- cbind(X, B_res$X_ffvd)
    Z <- cbind(Z, B_res$Z_ffvd)
    G <- c(G, B_res$G_ffvd)
  }

  for (column_index in rev(non_special_indices)) {
    X <- cbind(evals[[column_index]], X)
  }

  if (is.null(Z)) {
    stop("please use the 'lm' function from the 'stats' package", call. = TRUE)
  }

  X <- as.matrix(cbind(rep(1, lenght = nrow(X)), X))

  G <- add_zeros_to_G(G, nf, ncol(Z))

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

  intercept <- fit$b.fixed[1]

  Vp_tmp <- calculate_Vp_tmp(fit, non_special_indices, nf, n_coefs_f)

  theta_aux <- calculate_theta_aux(fit, non_special_indices, nf, l.f)

  covar_theta <- B_res$TMatrix %*% Vp_tmp %*% t(B_res$TMatrix)

  theta_ffvd <- B_res$TMatrix %*% theta_aux

  theta_no_functional <- calculate_theta_no_functional(fit, non_special_indices)

  theta_f <- calculate_theta_f(nf, l.f, fit, non_special_indices)
  Beta_ffvd <- calculate_beta_ffvd(nffvd, data, response, M, L_Phi, B_T, theta_ffvd, deglist)

  # ffvd_evals <- process_ffvd_evals(evals)
  ffvd_evals <- process_evals(evals, type = "ffvd")

  res <- list(
    fit = fit,
    Beta = Beta_ffvd,
    theta = theta_ffvd,
    covar_theta = covar_theta,
    M = M,
    ffvd_evals = ffvd_evals,
    theta_no_functional = theta_no_functional,
    theta_f = theta_f
  )

  class(res) <- "vd_fit"
  attr(res, "N") <- length(data[[response]])

  res
}



#' @export
summary.vd_fit <- function(object, ...) {
  base::summary(object$fit, ...)
}

#' Helper function to process ffvd evaluations
#'
#' @noRd
process_ffvd_evals <- function(evals) {
  ffvd_evals <- evals[grepl("ffvd", names(evals))]
  names <- paste0("ffvd_", gsub(".*?\\((.+?),.+", "\\1", names(ffvd_evals)))
  names(ffvd_evals) <- names
  ffvd_evals
}

#' Helper function to compute beta_ffvd
#'
#' @noRd
calculate_beta_ffvd <- function(nffvd, data, response, M, L_Phi, B_T, theta, deglist) {
  if (nffvd == 0) {
    return(NULL)
  }

  nffvd <- length(M)
  Beta_ffvd <- vector("list", nffvd)

  # Initialize Beta_ffvd matrices
  for (j in seq_len(nffvd)) {
    Beta_ffvd[[j]] <- matrix(
      NA,
      nrow = length(data[[response]]),
      ncol = max(M[[j]])
    )
  }

  # Calculate Beta_ffvd values
  theta_index <- 1
  for (j in seq_len(nffvd)) {
    for (i in seq_along(data[[response]])) {
      range_start <- M[[j]][i, 1]
      range_end <- M[[j]][i, 2]
      range <- range_start:range_end

      kron_product <- kronecker(
        L_Phi[[j]]$B[range, ],
        t(B_T[[j]]$B[i, ])
      )

      end_index <- theta_index + prod(deglist[[j]]) - 1
      Beta_ffvd[[j]][i, range] <- as.matrix(kron_product) %*% theta[theta_index:end_index]
    }
    theta_index <- theta_index + prod(deglist[[j]])
  }

  Beta_ffvd
}

#' Helper function to sum the number of column for all the matrices inside a list
#'
#' @noRd
sum_cols <- function(list, type) {
  sum(sapply(list, function(x) ncol(x$Xmat[[type]])))
}

#' Helper function to compute theta_aux
#'
#' @noRd
calculate_theta_aux <- function(fit, non_special_indices, nf, l.f = NULL) {
  # Determine case and calculate accordingly
  case_type <- list(
    has_non_special = length(non_special_indices) > 0,
    has_f = nf > 0
  )

  if (case_type$has_non_special && !case_type$has_f) {
    # Only non-special indices
    excluded_range <- 1:(length(non_special_indices) + 1)
    b_fixed_subset <- fit$b.fixed[-excluded_range]
    return(c(b_fixed_subset, fit$b.random))
  } else if (!case_type$has_non_special && case_type$has_f) {
    # Only f
    x_cols <- sum_cols(l.f, "X")
    z_cols <- sum_cols(l.f, "Z")

    b_fixed_subset <- fit$b.fixed[-c(1:x_cols)]
    b_random_subset <- fit$b.random[-c(1:z_cols)]
    return(c(b_fixed_subset, b_random_subset))
  } else if (case_type$has_non_special && case_type$has_f) {
    # Both f and non-special indices
    last_index_z <- sum_cols(l.f, "Z")
    x_cols <- sum_cols(l.f, "X")

    last_index_x <- seq(
      length(non_special_indices) + 2,
      by = 1,
      length.out = x_cols
    )

    b_fixed_subset <- fit$b.fixed[(max(last_index_x) + 1):length(fit$b.fixed)]
    b_random_subset <- fit$b.random[(last_index_z + 1):length(fit$b.random)]
    return(c(b_fixed_subset, b_random_subset))
  } else {
    # Neither f nor non-special indices
    return(c(fit$b.fixed[-1], fit$b.random))
  }
}

#' Helper function to compute theta_f
#'
#' @noRd
calculate_theta_f <- function(nf, l.f, fit, non_special_indices) {
  if (nf == 0) {
    return(NULL)
  }

  # Calculate total columns for X matrices and last index for random effects
  total_x_cols <- sum(sapply(l.f, function(x) ncol(x$Xmat$X)))
  last_random_index <- sum(sapply(l.f, function(x) ncol(x$Xmat$Z)))

  # Generate indices for fixed effects
  indices_b_fixed <- if (non_special_indices) {
    seq(
      length(non_special_indices) + 2,
      by = 1,
      length.out = total_x_cols
    )
  } else {
    1:total_x_cols
  }

  # Combine fixed and random effects
  c(
    fit$b.fixed[indices_b_fixed],
    fit$b.random[1:last_random_index]
  )
}

#' Helper function to add zeros to the different Gs
#'
#' @noRd
add_zeros_to_G <- function(G, nf, ncol_Z) {
  if (nf == 0) {
    return(G)
  }

  for (i in seq_along(G)) {
    side <- if (i <= nf) "right" else "left"
    G[[i]] <- add_zeros_to_side(G[[i]], ncol_Z, side = side)
  }

  G
}

#' Helper function to calculate theta_no_functional
#'
#' @noRd
calculate_theta_no_functional <- function(fit, non_special_indices) {
  if (length(non_special_indices) == 0) {
    return(NULL)
  }

  fit$b.fixed[2:(length(non_special_indices) + 1)]
}

#' Helper function to calculate Vp_tmp
#'
#' @noRd
calculate_Vp_tmp <- function(fit, non_special_indices, nf, n_coefs_f) {
  # Calculate total offset
  base_offset <- 1 # We always exclude at least the first row/column
  special_offset <- length(non_special_indices)
  f_offset <- if (nf > 0) n_coefs_f else 0

  total_offset <- base_offset + special_offset + f_offset

  # If the total offset is the base offset we only exclude the first row/column
  if (total_offset == 1) {
    return(fit$Vp[-1, -1])
  }

  exclude_indices <- 1:total_offset

  fit$Vp[-exclude_indices, -exclude_indices]
}
