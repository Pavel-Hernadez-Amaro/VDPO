library(refund)
library(fda)
library(mgcv)
library(SOP)

options(error=NULL)

########### Here we generate the data and set some fixed parameters

N=10 # NUMBER OF SUBJECTS

J=100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

R=1 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY

case=1 # THIS IS THE TOTAL NUMBER OF SCENARIOS SIMULATED


for (iter_out in 1:case) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE

  print(c("case = ",iter_out))

  for (iter in 1:R) {

    print(c("iter = ",iter))

    set.seed(1000+iter)

    M_1 = round(runif(N,1,J/2),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS

    M_2 = round(runif(N,M_1,J),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
    # M = rnegbin(N,24,1) # Para 1000 poner (240,2)

    M=cbind(M_1,M_2)

    # if (max(M)>J) {
    #   # print("si")
    #   M[which(M>J)]=J
    # }
    #
    #
    # if (min(M)<=10) {
    #   # print("si")
    #   # M[which(M<=10)]=round(runif(length(which(M<=10)),31,max(M)))
    #   M[which(M<=10)]=10
    # }

    M_diff=M[,2]-M[,1]+1

    T=max(M)

    t=min(M):max(M)

    # M=sort(M) # WE SORT THE DATA WITHOUT LOSS OF GENERALITY

    ############ HERE WE GENERATE THE FUNCTIONAL DATA

    X_se=matrix(NA,N,T) # NOISY

    X_s=matrix(NA,N,T) # NOT NOISY

    for (i in 1:N) {

      u=rnorm(1)

      temp=matrix(NA,10,T)

      for (k in 1:10) {

        v_i1=rnorm(1,0,4/k^2)
        v_i2=rnorm(1,0,4/k^2)

        temp[k,(M[i,1]:M[i,2])]=v_i1*sin(2*pi*k*(M[i,1]:M[i,2])/100)+v_i2*cos(2*pi*k*(M[i,1]:M[i,2])/100)
      }

      B=colSums(temp)[which(!is.na(colSums(temp)))]
      B=B+u

      X_s[i,M[i,1]:M[i,2]]=B
      X_se[i,M[i,1]:M[i,2]]=(B)+rnorm(M_diff[i],0,1) # WE ADD NOISE

    }

    # SOME SAVE CHECKS FOR UNWANTED NAs
    for (i in 1:T) {
      if (length(which(is.na(X_s[,i])))==N) {
        print(c("iteracion",l,"columna",i))
        stop("Hay columnas con NA")
      }}


    Beta=array(dim = c(N,T,4))
    nu=y=rep(0,N)

    for (i in 1:N) {

      # TRUE FUNCTIONAL COEFFICIENTS

      Beta[i,(M[i,1]:M[i,2]),1]=((10*t[(M[i,1]:M[i,2])]/M_diff[i])-5)/10
      Beta[i,(M[i,1]:M[i,2]),2]=((1-(2*M_diff[i]/T))*(5-40*((t[(M[i,1]:M[i,2])]/M_diff[i])-0.5)^2))/10
      Beta[i,(M[i,1]:M[i,2]),3]=(5-10*((M_diff[i]-t[(M[i,1]:M[i,2])])/T))/10
      Beta[i,(M[i,1]:M[i,2]),4]=(sin(2*pi*M_diff[i]/T)*(5-10*((M_diff[i]-t[(M[i,1]:M[i,2])])/T)))/10
      #


      # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)

      if (iter_out==8) {
        nu[i]=sum(X_se[i,]*Beta[i,,4],na.rm = 1)/(M_diff[i]) # NOISY
      }
      if (iter_out<=4) {
        nu[i]=sum(X_s[i,]*Beta[i,,iter_out],na.rm = 1)/(M_diff[i]) #NOT NOISY
      }
      if(iter_out>4 & iter_out<8){
        nu[i]=sum(X_se[i,]*Beta[i,,iter_out%%4],na.rm = 1)/M_diff[i] # NOISY
      }

    }


    y=nu+rnorm(N,sd = 1) # ADDING NOISE TO THE GAUSSIAN MODEL

    # y=rpois(N,exp(nu)) # POISSON MODEL ##### CAHNGE THE ERROR EQUATIONS IN LINE 477

    #     for (i in 1:N) {
    #
    #       # TRUE FUNCTIONAL COEFFICIENTS
    #
    #
    #
    #       # NUMBER OF BASIS FOR EVERY MARGINAL FUNCTION IN OUR APPROACH
    #       c1=25
    #       c2=25
    #       c3=25
    #
    #
    #       BB=Data2B_simpson(X_se, M, nbasis=c(c1,c2,c3),sub = 500, lim =c(min(M),max(M))) # HERE WE TRASFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL
    #       E=B2XZG(BB$B,c=c(c2,c3)) # HERE THE FUNCTION B2XZG() GENERATE THE NECESSARY MATRICES FOR THE MIXED MODEL
    #       res=XZG2theta(E$X, E$Z, E$G, E$T, y,family = gaussian()) # HERE THE MODEL COEFICIENT ARE ESTIMATED.
    #
    #
    #
    # }
  }}

data = data.frame(y = y)
data[["X"]] <- X_se
data[["Y"]] <- X_s

formula = y~ffvd(X)+ffvd(Y)

res <- VDFO(formula, data = data)
#
# plot(X_se, type = "l")
