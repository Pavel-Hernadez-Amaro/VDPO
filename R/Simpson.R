#' Simpson numeric integration
#'
#' Simpson numeric integration proposed in \href{http://dx.doi.org/10.9734/BJMCS/2016/23048}{Burden and Fairies (2016)}
#'
#' @param fdobj1 First functional object.
#' @param fdobj2 Second functional object.
#' @param fdobj3 Third functional object.
#' @param Lfdobj1 .
#' @param Lfdobj2 .
#' @param rng .
#' @param sub Number of intervals for the Simpson numerical integration method
#' divided by 2.
#' @param wtfd .
#'
#' @return .
#'
#' @noRd
Simpson <- function(fdobj1, fdobj2=NULL, fdobj3=NULL, Lfdobj1=fda::int2Lfd(0), Lfdobj2=fda::int2Lfd(0), rng, sub = 25, wtfd = 0) {

  #  Check FDOBJ1 and get no. replications and basis object

  result1   <- fdchk(fdobj1)
  nrep1     <- result1[[1]]
  fdobj1    <- result1[[2]]
  coef1     <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1     <- basisobj1$type
  range1    <- basisobj1$rangeval

  #  Default FDOBJ2 to a constant function, using a basis that matches
  #  that of FDOBJ1 if possible.

  if (is.null(fdobj2)) {
    tempfd    <- fdobj1
    tempbasis <- tempfd$basis
    temptype  <- tempbasis$type
    temprng   <- tempbasis$rangeval
    if (temptype == "bspline") {
      basis2 <- fda::create.bspline.basis(temprng, 1, 1)
    } else {
      if (temptype == "fourier") basis2 <- fda::create.fourier.basis(temprng, 1)
      else                       basis2 <- fda::create.constant.basis(temprng)
    }
    fdobj2 <- fda::fd(1,basis2)
  }

  #  Check FDOBJ2 and get no. replications and basis object

  result2   <- fdchk(fdobj2)
  nrep2     <- result2[[1]]
  fdobj2    <- result2[[2]]
  coef2     <- fdobj2$coefs
  basisobj2 <- fdobj2$basis
  type2     <- basisobj2$type
  range2    <- basisobj2$rangeval


  if (is.null(fdobj3)) {
    return(inprod(fdobj1, fdobj2, Lfdobj1=fda::int2Lfd(0), Lfdobj2=fda::int2Lfd(0),
                  rng = rng, wtfd = 0))
  }

  # check ranges

  if (rng[1] < range1[1] || rng[2] > range1[2]) stop(
    "Limits of integration are inadmissible.")

  #  check LFDOBJ1 and LFDOBJ2

  Lfdobj1 <- fda::int2Lfd(Lfdobj1)
  Lfdobj2 <- fda::int2Lfd(Lfdobj2)
  # Lfdobj3 <- int2Lfd(Lfdobj3)

  #  set iter

  iter <- 0

  # The default case, no multiplicities.

  rngvec <- rng

  #  check for any knot multiplicities in either argument

  knotmult <- numeric(0)

  if (type1 == "bspline") knotmult <- knotmultchk(basisobj1, knotmult)
  if (type2 == "bspline") knotmult <- knotmultchk(basisobj2, knotmult)

  #  Modify RNGVEC defining subinvervals if there are any
  #  knot multiplicities.

  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < rng[2]]
    rngvec   <- c(rng[1], knotmult, rng[2])
  }

  #  check for either coefficient array being zero

  if ((all(c(coef1) == 0) || all(c(coef2) == 0)))
    return(matrix(0,nrep1,nrep2*length(fdobj3)))

  # us: we need to define the number of intervals for the Simpson method
  n <- 2 * sub

  nrng <- length(rngvec)
  for (irng  in  2:nrng) {
    rngi <- c(rngvec[irng-1],rngvec[irng])

    #  change range so as to avoid being exactly on multiple knot values

    if (irng > 2   ) rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) rngi[2] <- rngi[2] - 1e-10

    width <- (rngi[2] - rngi[1])/n

    # the first iteration uses just the endpoints
    fx1 <- fda::eval.fd(rngi, fdobj1, Lfdobj1)
    fx2 <- fda::eval.fd(rngi, fdobj2, Lfdobj2)
    fx3 <- matrix(fdobj3,nrow = dim(fx2)[1],ncol = length(fdobj3),byrow = TRUE)
    fx2 <- Rten2(fx2,fx3)

    XI0 <- matrix(crossprod(fx1,fx2),nrep1,nrep2*length(fdobj3))

    XI1 <- 0
    XI2 <- 0

    for (i in 1:(n-1)) {

      ### THESE LINES OF CODE ARE THE CORE OF THE INNER PRODUCT. NOTICE HOW WE PERFORM THE INTEGRATION ONLY IN THE t VARIABLE BUT THEN WE CREATE THE BIDEMENSIONAL BASIS IN EVERY ITERATION

      x <- rngi[1] + i * width
      fx1 <- fda::eval.fd(x, fdobj1, Lfdobj1)
      fx2 <- fda::eval.fd(x, fdobj2, Lfdobj2)
      fx3 <- matrix(fdobj3,nrow = dim(fx2)[1],ncol = length(fdobj3),byrow = TRUE)
      fx2 <- Rten2(fx2,fx3)

      Fx <- matrix(crossprod(fx1,fx2),nrep1,nrep2*length(fdobj3))

      if (i %% 2 == 0) {
        XI2 <- XI2 + Fx
      } else {
        XI1 <- XI1 + Fx
      }

    }

    XI <- width * (XI0 + 2 * XI2 + 4 * XI1) / 3
  }

  if(length(dim(XI) == 2)) {
    #  coerce inprodmat to be nonsparse
    return(as.matrix(XI))
  } else {
    #  allow inprodmat to be sparse if it already is
    return(XI)
  }

}
