VDFO <- function(formula, data, family = stats::gaussian(), offset = NULL) {
  if (inherits(data, what = "data.frame")) {
    data <- as.data.frame(data)
  } else {
    stop("The data specified in the 'data' argument should be a data frame",
         call. = FALSE)
  }

  nobs <- nrow(data)
  if (is.null(offset))
    offset <- rep(0L, nobs)

  tf <- terms.formula(formula, specials = c("ffvd"))

  terms    <- attr(tf, "term.labels")
  nterms   <- length(terms)
  specials <- attr(tf, "specials") # indices for the special terms
  # specials <- any(unlist(lapply(specials_indices,
  #                        function (x) length(x) > 0)))


  if (attr(tf, "response")) {
    response <- as.character(attr(tf, "variables"))[2]
    variables <- as.character(attr(tf, "variables"))[-c(1,2)]

    # If the response exists, we need to lower in the index for every
    # special term

    specials <- vapply(specials,
                               function(x) x - 1,
                               FUN.VALUE = numeric(1))
  } else {
    variables <- as.character(attr(tf, "variables"))[-1]

    # No need to lower the index of the specials terms if the response exists
    specials <- vapply(specials,
                               function(x) if (length(x) > 0) x,
                               FUN.VALUE = numeric(1))
  }

  non_funcional_variables <- variables[-specials]
  functional_variables    <- variables[specials]

  # terms(formula)
  vdfoenv <- environment(formula)

  for (var in names(data))
    vdfoenv[[var]] <- data[[var]]

  vdfons  <- loadNamespace("VDPO")

  ffvd_output <- eval(parse(text = terms[1]), envir = vdfons, enclos = vdfoenv)

  fit <- sop.fit(
    X = ffvd_output[["X_ffvd"]], Z = ffvd_output[["Z_ffvd"]], G = ffvd_output[["G_ffvd"]],
    y = data[[response]], family = family,
    control = list(trace = FALSE), offset = offset
  )

  theta_aux <- c(fit$b.fixed,fit$b.random)
  theta <- ffvd_output[["TMatrix"]] %*% theta_aux

  covar_theta <- ffvd_output[["TMatrix"]] %*% fit$Vp %*% t(ffvd_output[["TMatrix"]])
  std_error_theta <- sqrt(diag(ffvd_output[["TMatrix"]] %*% fit$Vp %*% t(ffvd_output[["TMatrix"]])))

  Beta_ffvd <- matrix(NA, nrow = length(data[[response]]),
                      ncol = max(ffvd_output$M)) ##

  for (i in 1:length(data[[response]])) {
    Beta_ffvd[i, 1:ffvd_output$M[i]] <- as.matrix(kronecker(ffvd_output$L_Phi[[i]]$B, t(ffvd_output$B_T$B[i,]))) %*% theta
  }

  res <- list(
    fit = fit,
    theta = theta,
    Beta_ffvd = Beta_ffvd,
    M_ffvd = ffvd_output$M
  )

  class(res) <- "VDFO"
  attr(res, "N") <- length(data[[response]])
  res

}

