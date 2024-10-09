#' Estimation of the generalized additive functional regression models for
#' variable domain functional data
#'
#' The \code{vd_fit} function fits generalized additive functional regression models
#' for variable domain functional data.
#'
#' @param formula a formula object with at least one \code{ffvd} term.
#' @param data a \code{list} object containing the response variable
#' and the covariates.
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
#' \dontrun{
#' # VARIABLE DOMAIN FUNCTIONAL DATA EXAMPLE
#'
#' # load the example data
#' data <- VDPO::VDPO_example_vd
#'
#' # define a formula object that determines the model behavior
#' # note that this formula only uses one 'ffvd' term and that
#' # the 'nbasis' parameter is not the default one
#' formula <- y ~ ffvd(X_se, nbasis = c(10, 10, 10))
#'
#' # fit the model with the data and the formula
#' res <- VDPO(formula = formula, data = data)
#'
#' # important parameters of the model can be accessed as follows
#' res$Beta_ffvd         # variable domain functional coefficient
#' res$fit$fitted.values # estimated response variable
#'
#' # ------------------------------------------------------------------
#' # PARTIALLY OBSERVED FUNCTIONAL DATA EXAMPLE
#'
#' # load the example data
#' data <- VDPO::VDPO_example_po
#'
#' # define a formula object that determines the model behavior
#' # note that this formula only uses one 'ffpo' term
#' formula <- y ~ ffvd(X_se)
#'
#' # fit the model with the data and the formula
#' res <- VDPO(formula = formula, data = data)
#'
#' # important parameters of the model can be accessed as follows
#' res$theta_ffpo        # functional coefficient
#' res$fit$fitted.values # estimated response variable
#'}
#' @seealso \code{\link{ffvd}}, \code{\link{ffpo}}, \code{\link{ffpo_2d}}, \code{\link{add_grid}}
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
  nffvd    <- sum(grepl("\\bffvd\\(\\b", names(evals)))

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
    B_all   <- c()
    deglist <- vector(mode = "list", length = nffvd)
    L_Phi   <- vector(mode = "list", length = nffvd) # list of matrices
    B_T     <- vector(mode = "list", length = nffvd) # list of matrices
    M       <- vector(mode = "list", length = nffvd) # list of matrices
    # TMatrix <- vector(mode = "list", length = nffvd) # matrix

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


  # We need to add zeros to the right of the f Gs and to the left of the ffvd

  if (nf > 0) {
    for (i in seq_along(G)) {
      G[[i]] <- if (i <= nf) {
        add_zeros_to_side(G[[i]], ncol(Z), side = "right")
      } else {
        add_zeros_to_side(G[[i]], ncol(Z), side = "left")
      }
    }
  }
  # Adding zeros to the left for the ffvd
  # for (i in (nf+1):(2*(nf+nffvd)-1)) {
  #   G[[i]] <- c(rep(0, ncol(Z) - length(G[[i]])), G[[i]])
  # }

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

  Vp_tmp <- if (length(non_special_indices) > 0 && nf == 0) {
    # Only non special indices
    fit$Vp[-c(1:(length(non_special_indices) + 1)), -c(1:(length(non_special_indices) + 1))]

  } else if (length(non_special_indices) == 0 && nf > 0) {
    # Only f
    fit$Vp[-c(1:(n_coefs_f + 1)), -c(1:(n_coefs_f + 1))]

  } else if (length(non_special_indices) > 0 && nf > 0) {
    # f and non special indices
    fit$Vp[-c(1:(length(non_special_indices) + n_coefs_f + 1)), -c(1:(length(non_special_indices) + n_coefs_f + 1))]

  } else {
    fit$Vp[-1, -1]
  }

  theta_aux <- calculate_theta_aux(fit, non_special_indices, nf, l.f)

  covar_theta <- B_res$TMatrix %*% Vp_tmp %*% t(B_res$TMatrix)

  theta_ffvd <- B_res$TMatrix %*% theta_aux

  theta_no_functional <- if (length(non_special_indices)) {
    fit$b.fixed[2:(length(non_special_indices)+1)]
  } else {
    NULL
  }

  theta_f <- if (nf) {
    last_index <- sum(sapply(l.f, function(x) ncol(x$Xmat$Z)))

    indices_b_fixed <- if (non_special_indices) {
        seq(
          length(non_special_indices) + 2,
          by = 1,
          length.out = sum(sapply(l.f, function(x) ncol(x$Xmat$X)))
        )
    } else {
        1:sum(sapply(l.f, function(x) ncol(x$Xmat$X)))
    }

    c(fit$b.fixed[indices_b_fixed], fit$b.random[1:last_index])
  } else {
    NULL
  }

  Beta_ffvd <- if (nffvd > 0) {
    calculate_beta_ffvd(data, response, M, L_Phi, B_T, theta_ffvd, deglist)
  } else {
    NULL
  }

  ffvd_evals <- process_ffvd_evals(evals)

  res <- list(
    fit         = fit,
    Beta        = Beta_ffvd,
    theta       = theta_ffvd,
    covar_theta = covar_theta,
    M           = M,
    ffvd_evals  = ffvd_evals,
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
calculate_beta_ffvd <- function(data, response, M, L_Phi, B_T, theta, deglist) {
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

#' Helper function to sum columns
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
