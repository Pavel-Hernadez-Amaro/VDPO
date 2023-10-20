#' Inner product for fdobejct
#'
#' @param fdobj1 .
#' @param fdobj2 .
#' @param Lfdobj1 .
#' @param Lfdobj2 .
#' @param rng .
#' @param wtfd .
#'
#' @return .
#'
#' @noRd
inprod <- function (fdobj1, fdobj2 = NULL, Lfdobj1 = fda::int2Lfd(0), Lfdobj2 = fda::int2Lfd(0), rng = range1, wtfd = 0) {

  result1   <- fdchk(fdobj1)
  nrep1     <- result1[[1]]
  fdobj1    <- result1[[2]]
  coef1     <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1     <- basisobj1$type
  range1    <- basisobj1$rangeval

  if (is.null(fdobj2)) {
    tempfd <- fdobj1
    tempbasis <- tempfd$basis
    temptype <- tempbasis$type
    temprng <- tempbasis$rangeval
    if (temptype == "bspline") {
      basis2 <- fda::create.bspline.basis(temprng, 1, 1)
    }
    else {
      if (temptype == "fourier")
        basis2 <- fda::create.fourier.basis(temprng, 1)
      else basis2 <- fda::create.constant.basis(temprng)
    }
    fdobj2 <- fda::fd(1, basis2)
  }

  result2   <- fdchk(fdobj2)
  nrep2     <- result2[[1]]
  fdobj2    <- result2[[2]]
  coef2     <- fdobj2$coefs
  basisobj2 <- fdobj2$basis
  type2     <- basisobj2$type
  range2    <- basisobj2$rangeval

  if (rng[1] < range1[1] || rng[2] > range1[2])
    stop("Limits of integration are inadmissible.")

  # Relaxed version of the original version of the algorithm. Commented.

  # if (is.fd(fdobj1) && is.fd(fdobj2) && type1 == "bspline" &&
  #     type2 == "bspline" && is.eqbasis(basisobj1, basisobj2) &&
  #     is.integer(Lfdobj1) && is.integer(Lfdobj2) && length(basisobj1$dropind) ==
  #     0 && length(basisobj1$dropind) == 0 && wtfd == 0 && all(rng ==
  #                                                             range1)) {
  #   inprodmat <- inprod.bspline(fdobj1, fdobj2, Lfdobj1$nderiv,
  #                               Lfdobj2$nderiv)
  #   return(inprodmat)
  # }

  Lfdobj1 <- fda::int2Lfd(Lfdobj1)
  Lfdobj2 <- fda::int2Lfd(Lfdobj2)
  iter <- 0
  rngvec <- rng
  knotmult <- numeric(0)
  if (type1 == "bspline")
    knotmult <- knotmultchk(basisobj1, knotmult)
  if (type2 == "bspline")
    knotmult <- knotmultchk(basisobj2, knotmult)
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult <
                           rng[2]]
    rngvec <- c(rng[1], knotmult, rng[2])
  }
  if ((all(c(coef1) == 0) || all(c(coef2) == 0)))
    return(matrix(0, nrep1, nrep2))
  JMAX <- 15
  JMIN <- 5
  EPS <- 1e-04
  inprodmat <- matrix(0, nrep1, nrep2)
  nrng <- length(rngvec)
  for (irng in 2:nrng) {
    rngi <- c(rngvec[irng - 1], rngvec[irng])
    if (irng > 2)
      rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng)
      rngi[2] <- rngi[2] - 1e-10
    iter <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1, JMAXP)
    h[2] <- 0.25
    s <- array(0, c(JMAXP, nrep1, nrep2))
    sdim <- length(dim(s))
    fx1 <- fda::eval.fd(rngi, fdobj1, Lfdobj1)
    fx2 <- fda::eval.fd(rngi, fdobj2, Lfdobj2)

    if (!is.numeric(wtfd)) {
      wtd <- fda::eval.fd(rngi, wtfd, 0)
      fx2 <- matrix(wtd, dim(wtd)[1], dim(fx2)[2]) * fx2
    }
    s[1, , ] <- width * matrix(crossprod(fx1, fx2), nrep1, nrep2)/2
    tnm <- 0.5

    for (iter in 2:JMAX) {
      tnm <- tnm * 2
      if (iter == 2) {
        x <- mean(rngi)
      }
      else {
        del <- width/tnm
        x   <- seq(rngi[1] + del/2, rngi[2] - del/2, del)
      }

      fx1 <- fda::eval.fd(x, fdobj1, Lfdobj1)
      fx2 <- fda::eval.fd(x, fdobj2, Lfdobj2)
      if (!is.numeric(wtfd)) {
        wtd <- fda::eval.fd(wtfd, x, 0)
        fx2 <- matrix(wtd, dim(wtd)[1], dim(fx2)[2]) *
          fx2
      }
      chs <- width * matrix(crossprod(fx1, fx2), nrep1, nrep2)/tnm
      s[iter, , ] <- (s[iter - 1, , ] + chs)/2
      if (iter >= 5) {
        ind <- (iter - 4):iter
        ya <- s[ind, , ]
        ya <- array(ya, c(5, nrep1, nrep2))
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y <- ya[ns, , ]
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5 - m)) {
            ho <- xa[i]
            hp <- xa[i + m]
            w <- (cs[i + 1, , ] - ds[i, , ])/(ho - hp)
            ds[i, , ] <- hp * w
            cs[i, , ] <- ho * w
          }
          if (2 * ns < 5 - m) {
            dy <- cs[ns + 1, , ]
          }
          else {
            dy <- ds[ns, , ]
            ns <- ns - 1
          }
          y <- y + dy
        }
        ss <- y
        errval <- max(abs(dy))
        ssqval <- max(abs(ss))
        if (all(ssqval > 0)) {
          crit <- errval/ssqval
        }
        else {
          crit <- errval
        }
        if (crit < EPS && iter >= JMIN)
          break
      }
      s[iter + 1, , ] <- s[iter, , ]
      h[iter + 1] <- 0.25 * h[iter]
      if (iter == JMAX)
        warning("Failure to converge.")
    }
    inprodmat <- inprodmat + ss
  }
  if (length(dim(inprodmat) == 2)) {
    return(as.matrix(inprodmat))
  }
  else {
    return(inprodmat)
  }
}
