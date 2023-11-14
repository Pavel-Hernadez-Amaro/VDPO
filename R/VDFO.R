VDFO <- function(formula, data, family = stats::gaussian(), offset = NULL) {
  if (inherits(data, what = "data.frame")) {
    data <- as.data.frame(data)
  } else {
    stop("The data specified in the 'data' argument should be a data frame",
         call. = FALSE)
  }

  vdfoenv <- environment(formula)
  vdfons  <- loadNamespace("VDPO")

  for (var in names(data)) {
    vdfoenv[[var]] <- data[[var]]
  }

  nobs <- nrow(data)
  if (is.null(offset)) {
    offset <- rep(0L, nobs)
  }

  tf <- stats::terms.formula(formula, specials = c("ffvd"))

  terms    <- attr(tf, "term.labels")
  nterms   <- length(terms)
  specials_indices <- attr(tf, "specials") # indices for the special terms

  if (attr(tf, "response")) {
    response <- as.character(attr(tf, "variables"))[2]
    variables <- as.character(attr(tf, "variables"))[-c(1,2)]

    # If the response exists, we need to lower in the index for every
    # special term

    specials_indices <- lapply(specials_indices,
                       function(x) if (is.null(x)) NA else x - 1
                      )
  } else {
    variables <- as.character(attr(tf, "variables"))[-1]

    # No need to lower the index of the specials terms if the response exists
    specials_indices <- lapply(specials_indices,
                      function(x) if (is.null(x)) NA else x,
                      )
  }


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

  X_ffvd  <- c() # matrix
  Z_ffvd  <- c() # matrix
  G_ffvd  <- list() # list
  L_Phi   <- vector(mode = "list", length = sum(grepl("ffvd", names(evals)))) # list of matrices
  B_T     <- vector(mode = "list", length = sum(grepl("ffvd", names(evals)))) # matrix
  M       <- vector(mode = "list", length = sum(grepl("ffvd", names(evals)))) # matrix
  TMatrix <- vector(mode = "list", length = sum(grepl("ffvd", names(evals)))) # matrix

  lffvd_counter <- 1
  for (ffvd_evaluation in evals[grepl("ffvd", names(evals))]) {
    X_ffvd <- cbind(X_ffvd, ffvd_evaluation[["X_ffvd"]])
    Z_ffvd <- cbind(Z_ffvd, ffvd_evaluation[["Z_ffvd"]])
    G_ffvd <- append(G_ffvd, ffvd_evaluation[["G_ffvd"]])

    L_Phi[[lffvd_counter]]   <- ffvd_evaluation[["L_Phi"]]
    B_T[[lffvd_counter]]     <- ffvd_evaluation[["B_T"]]
    M[[lffvd_counter]]       <- ffvd_evaluation[["M"]]
    TMatrix[[lffvd_counter]] <- ffvd_evaluation[["TMatrix"]]

    lffvd_counter <- lffvd_counter + 1

  }

  fit <- sop.fit(
    X = X_ffvd, Z = Z_ffvd, G = G_ffvd,
    y = data[[response]], family = family,
    control = list(trace = FALSE), offset = offset
  )

  theta_aux <- c(fit$b.fixed, fit$b.random) ##
  theta <- lapply(TMatrix, function (x) x %*% theta_aux) # x is a TMatrix
  # theta <- evals[["TMatrix"]] %*% theta_aux

  # covar_theta <- evals[["TMatrix"]] %*% fit$Vp %*% t(evals[["TMatrix"]])
  covar_theta <- lapply(TMatrix, function(x) x %*% fit$Vp %*% t(x)) ##
  # std_error_theta <- sqrt(diag(evals[["TMatrix"]] %*% fit$Vp %*% t(evals[["TMatrix"]])))
  std_error_theta <- lapply(TMatrix, function(x) sqrt(diag(x %*% fit$Vp %*% t(x)))) ##

  # Beta_ffvd <- matrix(NA, nrow = length(data[[response]]),
  #                     ncol = max(evals$M)) ##
  Beta_ffvd <- lapply(M,
                      function (x) matrix(NA, nrow = length(data[[response]]),
                                          ncol = max(x))
                      )
  for (j in 1:sum(grepl("ffvd", names(evals)))) {
    for (i in 1:length(data[[response]])) {
      Beta_ffvd[[j]][i, M[[j]][i, 1]:M[[j]][i, 2]] <- as.matrix(kronecker(L_Phi[[j]][[i]]$B, t(B_T[[j]]$B[i,]))) %*% theta[[j]]
    }
  }
  res <- list(
    fit       = fit,
    theta     = theta,
    Beta_ffvd = Beta_ffvd,
    M_ffvd    = evals$M
  )

  class(res) <- "VDFO"
  attr(res, "N") <- length(data[[response]])
  res

}

