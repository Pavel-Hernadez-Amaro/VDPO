VDFO <- function(formula, data, family = stats::gaussian()) {
  if (inherits(data, what = "dataframe")) {
    data <- as.data.frame(data)
  } else {
    stop("The data specified in the 'data' argument should be a data frame",
         call. = FALSE)
  }

  tf <- terms.formula(formula, specials = c("ffvd"))

  terms    <- attr(tf, "term.labels")
  nterms   <- length(terms)
  specials <- attr(tf, "specials") # indices for the special terms

  if (attr(tf, "response")) {
    response <- as.character(attr(tf, "variables"))[2]
    variables <- as.character(attr(tf, "variables"))[-c(1,2)]

    # If the response exists, we need to lower in the index for every
    # special term

    specials_indices <- vapply(specials,
                               function(x) x - 1,
                               FUN.VALUE = numeric(1))
  } else {
    variables <- as.character(attr(tf, "variables"))[-1]

    # No need to lower the index of the specials terms if the response exists
    specials_indices <- vapply(specials,
                               function(x) x,
                               FUN.VALUE = numeric(1))
  }

  non_funcional_variables <- variables[-specials_indices]
  functional_variables    <- variables[specials_indices]

  terms(formula)
  vdfoenv <- environment(formula)
  vdfons  <- loadNamespace("VDPO")

  eval(parse(terms[1]), envir = vdfons, enclos = vdfoenv)





}

ffvd <- function(x) {
  ifelse(x<2, 10, 20)
}
