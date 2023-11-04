VDFO <- function(formula, data, family = stats::gaussian()) {
  if (inherits(data, what = "dataframe")) {
    data <- as.data.frame(data)
  } else {
    stop("The data specified in the 'data' argument should be a data frame",
         call. = FALSE)
  }

  tf <- terms.formula(formula, specials = c("ffvd"))

  terms <- attr(tf, "term.labels")
  nterms <- length(terms)
  specials <- attr(tf, "specials") # indices for the special terms

  if (attr(tf, "response")) {
    response <- as.character(attr(tf, "variables"))[2]
    variables <- as.character(attr(tf, "variables"))[-c(1,2)]

    specials$ffvd <- specials$ffvd - 1

    specials_indices <- c(specials$ffvd)
  } else {
    variables <- as.character(attr(tf, "variables"))[-1]
  }

  nfvariables <- variables[-specials_indices]
  fvariables <- variables[specials_indices]

  fvariables_names <- lapply(fvariables,
                             function (x) strsplit(
                               strsplit(x, "\\(")[[1]][[2]],
                               "\\)"
                               )[[1]][[1]]
  )

  fvariables_functions <- lapply(fvariables,
                                 function (x) strsplit(x, "\\(")[[1]][[1]])



  # ALL THE COMPUTATION




}

ffvd <- function(x) {
  ifelse(x<2, 10, 20)
}

f()

formula <- x ~ a + b + ffvd(c)
formula <- ~ a + b + ffvd(c)


sop()

m1 <- matrix(1:10,ncol=2)
m2 <- matrix(5:14,ncol=2)
dd <- list(m1,m2)

m <- matrix(7:10, nrow = 2)

data = data.frame(
  a = 1:2,
  b = 3:4,
  y = 5:6
)

data[["c"]] <- m

