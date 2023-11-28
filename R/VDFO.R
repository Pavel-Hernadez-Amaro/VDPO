VDFO <- function(formula, data, family = stats::gaussian(), offset = NULL) {
  if (inherits(data, what = "data.frame")) {
    data <- as.data.frame(data)
  } else {
    stop("The data specified in the 'data' argument should be a data frame",
         call. = FALSE)
  }

  vdfoenv <- environment(formula)
  vdfons  <- loadNamespace("VDPO")

  for (var in names(data))
    vdfoenv[[var]] <- data[[var]]

  nobs <- nrow(data)

  if (is.null(offset))
    offset <- rep(0L, nobs)

  tf <- stats::terms.formula(formula, specials = c("ffvd"))

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
  non_special_indices <- setdiff(all_indices, specials_indices)

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
  nffvd <- sum(grepl("ffvd", names(evals)))

  if (nffvd > 0) {
    B_all   <- c()
    deglist <- vector(mode = "list", length = nffvd)

    X_ffvd <- c() # matrix
    Z_ffvd <- c() # matrix
    G_ffvd <- list() # list
    L_Phi  <- vector(mode = "list", length = nffvd) # list of matrices
    B_T    <- vector(mode = "list", length = nffvd) # list of matrices
    M      <- vector(mode = "list", length = nffvd) # list of matrices
    # TMatrix <- vector(mode = "list", length = nffvd) # matrix

    ffvd_counter <- 0
    for (ffvd_evaluation in evals[grepl("ffvd", names(evals))]) {
      ffvd_counter <- ffvd_counter + 1

      B_all <- cbind(B_all, ffvd_evaluation[["B_ffvd"]])
      deglist[[ffvd_counter]] <- ffvd_evaluation[["nbasis"]][2:3]

      # X_ffvd <- cbind(X_ffvd, ffvd_evaluation[["X_ffvd"]])
      # Z_ffvd <- cbind(Z_ffvd, ffvd_evaluation[["Z_ffvd"]])
      # # G_ffvd <- append(G_ffvd, ffvd_evaluation[["G_ffvd"]])
      # G_ffvd <- append(G_ffvd,
      #                  list(c(ffvd_evaluation[["G_ffvd"]][[1]], rep(0, length(ffvd_evaluation[["G_ffvd"]][[2]]))),
      #                       c(rep(0, length(ffvd_evaluation[["G_ffvd"]][[1]])), ffvd_evaluation[["G_ffvd"]][[2]]))
      #                  )
      #
      L_Phi[[ffvd_counter]] <- ffvd_evaluation[["L_Phi"]]
      B_T[[ffvd_counter]]   <- ffvd_evaluation[["B_T"]]
      M[[ffvd_counter]]     <- ffvd_evaluation[["M"]]
      # TMatrix[[ffvd_counter]] <- ffvd_evaluation[["TMatrix"]]
    }
  }

  foo <- B2XZG(B_all, deglist)

  for (column_index in rev(non_special_indices))
    foo$X_ffvd <- cbind(evals[[column_index]], foo$X_ffvd)

  fit <- sop.fit(
    X = foo$X_ffvd, Z = foo$Z_ffvd, G = foo$G_ffvd,
    y = data[[response]], family = family,
    control = list(trace = FALSE), offset = offset
  )


  if (nrow(fit$Vp) == nrow(foo$TMatrix)) {
    theta_aux       <- c(fit$b.fixed,fit$b.random)
    covar_theta     <- foo$TMatrix %*% fit$Vp %*% t(foo$TMatrix)
    std_error_theta <- sqrt(diag(foo$TMatrix %*% fit$Vp %*% t(foo$TMatrix)))
  } else {
    aux             <- nrow(fit$Vp) - nrow(foo$TMatrix)
    theta_aux       <- c(fit$b.fixed[-(1:aux)], fit$b.random)
    covar_theta     <- foo$TMatrix %*% fit$Vp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)
    std_error_theta <- sqrt(diag(foo$TMatrix %*% fit$Vp[-(1:aux), -(1:aux)] %*% t(foo$TMatrix)))
    std_error_nf    <- sqrt(diag(fit$Vp[1:aux, 1:aux]))
    WALD            <-  (fit$b.fixed[(1:aux)] / std_error_nf)
    p_values        <- 2 * pnorm(abs(WALD), lower.tail = FALSE)
  }

  theta <- foo$TMatrix %*% theta_aux

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

  if (nrow(fit$Vp) == nrow(foo$TMatrix)) {
    res <- list(
      fit              = fit,
      theta            = theta,
      covar_theta      = covar_theta,
      std_error_theta  = std_error_theta,
      Beta_ffvd        = Beta_ffvd,
      M_ffvd           = M
    )
  } else { # if there are nf variables
    res <- list(
      fit             = fit,
      theta_nf        = fit$b.fixed[(1:aux)],
      theta           = theta,
      covar_theta     = covar_theta,
      std_error_theta = std_error_theta,
      std_error_nf    = std_error_nf,
      p_values        = p_values,
      Beta_ffvd       = Beta_ffvd,
      M_ffvd          = M
    )
  }



  class(res) <- "VDFO"
  attr(res, "N") <- length(data[[response]])
  res
}

