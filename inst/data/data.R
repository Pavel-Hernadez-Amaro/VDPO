response_integrator <- function(..., fx, fb, sub = 50, index = 1) {
  if (!requireNamespace("matrixcalc"))
    stop("package 'matrixcalc' is required for this functionality", call. = FALSE)

  dots <- list(...)
  N <- dots[["N"]]
  M <- dots[["M"]]

  n_y <- 2 * sub
  x <- seq(min(M), max(M), by = n_y)

  dg <- fx(N, M, precision = n_y)
  bg <- fb(N, M, precision = n_y)

  h <- seq(min(M), max(M), by = ncol(dg$X_se))  # if y_b and y_a are functions of x_observations

  simp_w_y <- rep(1,ncol(dg$X_se))
  even_y <- seq(2,ncol(dg$X_se)-1,2)
  odd_y <- seq(3,ncol(dg$X_se)-1,2)

  simp_w_y[even_y] <- 2
  simp_w_y[odd_y] <- 4

  Sim_w_x_y <- (h / 3) * simp_w_y

  dg$X_se <- ifelse(is.na(dg$X_se),
                    0, dg$X_se)
  bg[, , index] <- ifelse(is.na(bg[, , index]),
                          0, bg[, , index])


  sum(dg$X_se %*% diag(Sim_w_x_y) %*% t(bg[, , index]))

  # as.double(t(matrixcalc::vec(dg$X_se)) %*% diag(Sim_w_x_y) %*% matrixcalc::vec(bg[, , index]))

}


data_generator <- function(N, M, precision = 1) {
  precision <- 1 / precision
  sj <- seq(min(M), max(M), by = precision)
  X_se=matrix(NA,N,length(sj)) # NOISY
  X_s=matrix(NA,N,length(sj)) # NOT NOISY

  for (i in 1:N) {
    s <- seq(M[i,1], M[i,2], by = precision)
    u <- rnorm(1)

    temp <- matrix(NA, 10, length(s))


    for (k in 1:10) {

      v_i1=rnorm(1,0,4/k^2)
      v_i2=rnorm(1,0,4/k^2)

      temp[k, ] <-
        v_i1*sin(2*pi*k* s / 100)+v_i2*cos(2*pi*k*s/100)
    }

    B=colSums(temp)[which(!is.na(colSums(temp)))]
    B=B+u

    X_s[i,s]=B
    X_se[i,s]=(B)+rnorm(length(B), 0, 1) # WE ADD NOISE

  }

  list(
    X_s  = X_s,
    X_se = X_se
  )
}

beta_generator <- function(N, M, precision = 1) {
  precision <- 1 / precision
  rangeM <- seq(min(M), max(M), by = precision)
  maxM <- max(M)
  Beta <- array(dim = c(N,length(rangeM),4))
  nu=y=rep(0,N)

  for (i in 1:N) {
    s <- seq(M[i,1], M[i,2], by = precision)

    # TRUE FUNCTIONAL COEFFICIENTS

    Beta[i, s, 1]=((10*rangeM[s]/precision)-5)/10
    Beta[i, s, 2]=((1-(2*precision/maxM))*(5-40*((rangeM[s]/precision)-0.5)^2))/10
    Beta[i, s, 3]=(5-10*((precision-rangeM[s])/maxM))/10
    Beta[i, s, 4]=(sin(2*pi*precision/maxM)*(5-10*((precision-rangeM[s])/maxM)))/10

  }

  Beta
}
