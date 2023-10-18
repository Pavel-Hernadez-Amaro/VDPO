#' This Rten2 function, is, atm, a mystery to me.
#'
#' Auxiliary function found in some functions of the \href{https://cran.r-project.org/web/packages/sommer/index.html}{sommer} package.
#' All credits to them.
#'
#' @param X1 .
#' @param X2 .
#'
#' @references All credits to the \href{https://cran.r-project.org/web/packages/sommer/index.html}{sommer} package authors.
#'
Rten2 <- function(X1,X2) {
  one.1 <- matrix(1,1,ncol(X1))
  one.2 <- matrix(1,1,ncol(X2))
  kronecker(X1,one.2)*kronecker(one.1,X2)
}




#' @references All credits to the \href{https://cran.r-project.org/web/packages/fda/index.html}{fda} package authors.
fdchk <- function (fdobj) {
  if (inherits(fdobj, "fd")) {
    coef <- fdobj$coefs
  }
  else {
    if (inherits(fdobj, "basisfd")) {
      coef <- diag(rep(1, fdobj$nbasis - length(fdobj$dropind)))
      fdobj <- fda::fd(coef, fdobj)
    }
    else {
      stop("FDOBJ is not an FD object.")
    }
  }
  coefd <- dim(as.matrix(coef))
  if (length(coefd) > 2)
    stop("Functional data object must be univariate")
  nrep <- coefd[2]
  basisobj <- fdobj$basis
  return(list(nrep, fdobj))
}

#' @references All credits to the \href{https://cran.r-project.org/web/packages/fda/index.html}{fda} package authors.
knotmultchk <- function (basisobj, knotmult) {
  type <- basisobj$type
  if (type == "bspline") {
    params <- basisobj$params
    nparams <- length(params)
    norder <- basisobj$nbasis - nparams
    if (norder == 1) {
      knotmult <- c(knotmult, params)
    }
    else {
      if (nparams > 1) {
        for (i in 2:nparams) if (params[i] == params[i -
                                                     1])
          knotmult <- c(knotmult, params[i])
      }
    }
  }
  return(knotmult)
}
