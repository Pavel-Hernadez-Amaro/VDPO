VDFO <- function(formula, data, family = stats::gaussian(), offset = NULL) {
  if (inherits(formula, "character")) {
    formula <- stats::as.formula(formula)
  }
  if (inherits(data, what = "data.frame")) {
    data <- as.data.frame(data)
  } else {
    stop("The data specified in the 'data' argument should be a data frame",
         call. = FALSE)
  }

  na_indices <- c()
  for (name in names(data)) {
    if (is.null(dim(data[[name]]))) {
      na_indices <- c(na_indices, which(is.na(data[[name]])))
    }
  }

  na_indices <- unique(na_indices)

  if (length(na_indices) != 0)
    data <- data[-na_indices, ]

  vdfoenv <- environment(formula)
  vdfons  <- loadNamespace("VDPO")

  for (var in names(data))
    vdfoenv[[var]] <- data[[var]]

  nobs <- nrow(data)

  if (is.null(offset))
    offset <- rep(0L, nobs)

  tf <- stats::terms.formula(formula, specials = c("ffvd", "f"))

  terms  <- attr(tf, "term.labels")
  nterms <- length(terms)
  specials_indices <- attr(tf, "specials") # indices for the special terms

  if (attr(tf, "response")) {
    response  <- as.character(attr(tf, "variables"))[2]
    variables <- as.character(attr(tf, "variables"))[-c(1,2)]

    # If the response exists, we need to lower in the index for every
    # special term

    specials_indices <- lapply(specials_indices,
                               function(x) if (is.null(x)) NA else x - 1)
  } else {
    variables <- as.character(attr(tf, "variables"))[-1]

    # No need to lower the index of the specials terms if the response exists
    specials_indices <- lapply(specials_indices,
                               function(x) if (is.null(x)) NA else x)
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
    function(term) eval(parse(text = term), envir = vdfons, enclos = vdfoenv)
  )
  names(evals) <- terms
  nffvd <- sum(grepl("\\bffvd\\(\\b", names(evals)))
  ## TODO add the f function to the package namespace

  if (nffvd == 0) {
    stop("This function should be used with at least one 'ffvd' term.")
  }

  nf <- sum(grepl("\\bf\\(\\b", names(evals)))

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
      aa <- switch(as.character(l.f[[f_index]]$dim), `3` = {
        l.f[[f_index]]$Xmat <- SOP:::construct.3D.pspline(form,
                                                          data)
      }, `2` = {
        l.f[[f_index]]$Xmat <- SOP:::construct.2D.pspline(form,
                                                          data)
      }, `1` = {
        l.f[[f_index]]$Xmat <- SOP:::construct.1D.pspline(form,
                                                          data)
      })
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
      G <- SOP:::construct.capital.lambda(G)
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
    for (ffvd_evaluation in evals[grepl("ffvd", names(evals))]) {
      ffvd_counter <- ffvd_counter + 1

      B_all <- cbind(B_all, ffvd_evaluation[["B_ffvd"]])
      deglist[[ffvd_counter]] <- ffvd_evaluation[["nbasis"]][2:3]

      L_Phi[[ffvd_counter]] <- ffvd_evaluation[["L_Phi"]]
      B_T[[ffvd_counter]]   <- ffvd_evaluation[["B_T"]]
      M[[ffvd_counter]]     <- ffvd_evaluation[["M"]]
      # TMatrix[[ffvd_counter]] <- ffvd_evaluation[["TMatrix"]]
    }
    foo <- B2XZG(B_all, deglist)

    X <- cbind(X, foo$X_ffvd)
    Z <- cbind(Z, foo$Z_ffvd)
    G <- c(G, foo$G_ffvd)
  }



  for (column_index in rev(non_special_indices))
    X <- cbind(evals[[column_index]], X)

  if (is.null(Z))
    stop("check the 'lm' function from the 'stats' package.", call. = TRUE)

  X <- as.matrix(cbind(rep(1, lenght = nrow(X)), X))


  # We need to add zeros to the right of the f Gs and to to the left of the ffvd

  if (nf > 0) {
    for (i in 1:nf) {
      G[[i]] <- c(G[[i]], rep(0, ncol(Z) - length(G[[i]])))
    }
  }
  # Adding zeros to the right for the ffvd
  for (i in (nf+1):(2*(nf+nffvd)-1)) {
    G[[i]] <- c(rep(0, ncol(Z) - length(G[[i]])), G[[i]])
  }


  fit <- sop.fit(
    X = X, Z = Z, G = G,
    y = data[[response]], family = family,
    control = list(trace = FALSE), offset = offset
  )

  intercept <- fit$b.fixed[1]
  b_fixed_tmp <- fit$b.fixed[-1]
  Vp_tmp <- fit$Vp[-1,-1]

  if (length(non_special_indices) == 0 && nf == 0) {
    # only ffvd

    theta_ffvd      <- c( b_fixed_tmp,fit$b.random)
    covar_theta     <- foo$TMatrix %*%  Vp_tmp %*% t(foo$TMatrix)
    std_error_theta <- sqrt(diag(foo$TMatrix %*%  Vp_tmp %*% t(foo$TMatrix)))
  } else if (length(non_special_indices) > 0 && nf == 0) {
    # ffvd and x
    aux             <- nrow( Vp_tmp) - nrow(foo$TMatrix)
    theta_ffvd      <- c( b_fixed_tmp[-(1:aux)], fit$b.random)
    covar_theta     <- foo$TMatrix %*%  Vp_tmp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)
    std_error_theta <- sqrt(diag(foo$TMatrix %*%  Vp_tmp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)))
    std_error_nf    <- sqrt(diag( Vp_tmp[1:aux, 1:aux]))
    WALD            <-  ( b_fixed_tmp[(1:aux)] / std_error_nf)
    p_values        <- 2 * stats::pnorm(abs(WALD), lower.tail = FALSE)
  } else if (length(non_special_indices) != 0 && nf != 0) {
    # f, ffvd and x
    aux <- nrow( Vp_tmp) - nrow(foo$TMatrix)
    theta_f <- c( b_fixed_tmp[(aux + 1):(aux + aux_sum)], fit$b.random[1:aux_sum])
    theta_ffvd <- c( b_fixed_tmp[(aux + aux_sum + 1):length( b_fixed_tmp)], fit$b.random[(aux_sum + 1):length(fit$b.random)])
    # theta_aux here is theta_ffvd
    covar_theta     <- foo$TMatrix %*%  Vp_tmp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)
    std_error_theta <- sqrt(diag(foo$TMatrix %*%  Vp_tmp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)))

    std_error_nf    <- sqrt(diag( Vp_tmp[1:aux, 1:aux]))
    WALD            <-  ( b_fixed_tmp[(1:aux)] / std_error_nf)
    p_values        <- 2 * stats::pnorm(abs(WALD), lower.tail = FALSE)
  } else if (length(non_special_indices) == 0 && nf > 0) {
    # only ffvd and f

    aux <- nrow( Vp_tmp) - nrow(foo$TMatrix)
    theta_f <- c( fit$b.fixed[1], fit$b.random[1:aux])
    theta_ffvd <- c( fit$b.fixed[-1], fit$b.random[-(1:aux)])
    # theta_aux here is theta_ffvd
    covar_theta     <- foo$TMatrix %*%  Vp_tmp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)
    std_error_theta <- sqrt(diag(foo$TMatrix %*%  Vp_tmp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)))

  }


  theta <- foo$TMatrix %*% theta_ffvd

  Beta_ffvd <- lapply(M,
                      function (x) matrix(NA,
                                          nrow = length(data[[response]]),
                                          ncol = max(x)))
  it <- 1
  for (j in 1:nffvd) {
    for (i in 1:length(data[[response]])) {
      Beta_ffvd[[j]][i, M[[j]][i, 1]:M[[j]][i, 2]] <-
        as.matrix(kronecker(L_Phi[[j]][[i]]$B, t(B_T[[j]]$B[i,]))) %*% theta[it:(it + prod(deglist[[j]]) - 1)]
    }
    it <- it + prod(deglist[[j]])
  }

  if (length(non_special_indices) > 0 && nf > 0) {
    # ffvd and f and x

    res <- list(
      fit             = fit,
      theta_nf        = fit$b.fixed[(2:(aux+1))],
      theta_ffvd      = theta,
      theta_f         = theta_f,
      covar_theta     = covar_theta,
      std_error_theta = std_error_theta,
      std_error_nf    = std_error_nf,
      p_values        = p_values,
      Beta_ffvd       = Beta_ffvd,
      M_ffvd          = M
    )
  } else if (length(non_special_indices) > 0 && nf == 0) {
    # ffvd and x

    res <- list(
      fit             = fit,
      theta_nf        = fit$b.fixed[(2:(aux+1))],
      theta_ffvd      = theta,
      covar_theta     = covar_theta,
      std_error_theta = std_error_theta,
      std_error_nf    = std_error_nf,
      p_values        = p_values,
      Beta_ffvd       = Beta_ffvd,
      M_ffvd          = M
    )
  } else if(length(non_special_indices) == 0 && nf > 0) {
    # only ffvd and f

    res <- list(
      fit             = fit,
      theta_ffvd      = theta,
      theta_f         = theta_f,
      covar_theta     = covar_theta,
      std_error_theta = std_error_theta,
      Beta_ffvd       = Beta_ffvd,
      M_ffvd          = M
    )
  } else {
    # only ffvd

    res <- list(
      fit              = fit,
      theta_ffvd       = theta,
      covar_theta      = covar_theta,
      std_error_theta  = std_error_theta,
      Beta_ffvd        = Beta_ffvd,
      M_ffvd           = M
    )
  }

  class(res) <- "VDFO"
  attr(res, "N") <- length(data[[response]])
  res
}



#' @export
summary.VDFO <- function(object, ...) {
  base::summary(object$fit, ...)
}
