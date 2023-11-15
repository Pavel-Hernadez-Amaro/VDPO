# THIS FUNCTION GENERATE THE MATRICES X Z G TO BE USED IN THE SOP.FIT FOR THE 1D CASE WHEN WE HAVE TWO VARIABLES, SEE LEE AND DURBÁN REGUERA (2010)
B2XZG_1d_2_var_intercept <-function (B, pord=c(2,2),c=c(20,20)) {

  # with intercept

  c1=c[1]

  D_1 = diff(diag(c1), differences=pord[1])

  P1.svd = svd(crossprod(D_1))

  U_1s = (P1.svd$u)[,1:(c1-pord[1])] # eigenvectors
  U_1n=((P1.svd$u)[,-(1:(c1-pord[1]))])
  d1 = (P1.svd$d)[1:(c1-pord[1])]  # eigenvalues

  # T_n1=as.matrix((U_1n)[,-1])
  T_n1=as.matrix(U_1n)

  T_s1=U_1s

  c2=c[2]

  D_2 = diff(diag(c2), differences=pord[2])

  P2.svd = svd(crossprod(D_2))

  U_2s = (P2.svd$u)[,1:(c2-pord[2])] # eigenvectors
  U_2n=((P2.svd$u)[,-(1:(c2-pord[2]))])
  d2 = (P2.svd$d)[1:(c2-pord[2])]  # eigenvalues

  # T_n2=as.matrix((U_2n)[,-1])
  T_n2=as.matrix(U_2n)

  T_s2=U_2s

  # P=matrix(0,nrow = (1+dim(crossprod(D_1))[1]+dim(crossprod(D_2))[1]),ncol=(1+dim(crossprod(D_1))[1]+dim(crossprod(D_2))[1]))
  # P[2:(c1+1),2:(c1+1)]=crossprod(D_1)
  # P[(c1+2):(c1+c2+1),(c1+2):(c1+c2+1)]=crossprod(D_2)
  #
  # aux=svd(P)
  #
  # U_1=aux$u[2:(c1+1),2:(c1+1)]
  # U_1n=U_1[,-(1:(c1-pord[1]))]
  # U_1s=U_1[,1:(c1-pord[1])]
  #
  # U_2=aux$u[(c1+2):(c1+c2+1),(c1+2):(c1+c2+1)]
  # U_2n=U_2[,-(1:(c2-pord[2]))]
  # U_2s=U_2[,1:(c2-pord[2])]
  #
  # U_armada=matrix(0,nrow = dim(P)[1],ncol=dim(P)[2])
  # U_armada[2:(c1+1),2:(c1+1)]=P1.svd$u
  # U_armada[(c1+2):(c1+c2+1),(c1+2):(c1+c2+1)]=P2.svd$u
  #
  # all.equal(aux$u,U_armada)
  #
  #
  # all.equal(aux$d,c(P1.svd$d,P2.svd$d))



  T_n=rbind(cbind(T_n1,matrix(0,nrow = dim(T_n1)[1],ncol=dim(T_n2)[2])),cbind(matrix(0,nrow = dim(T_n2)[1],ncol=dim(T_n1)[2]),T_n2))
  T_n=rbind(diag(dim(T_n)[2]+1)[1,],cbind(0,T_n))

  T_s=rbind(cbind(T_s1,matrix(0,nrow = dim(T_s1)[1],ncol=dim(T_s2)[2])),cbind(matrix(0,nrow = dim(T_s2)[1],ncol=dim(T_s1)[2]),T_s2))
  T_s=rbind(0,T_s)

  # T=cbind(T_n,T_s)

  T=cbind(T_n[,-c(pord[1],pord[1]+pord[2])],T_s)  # THIS ONLY WORKS FOR pord=2 IN GENERAL IS THE SECOND COLUMN AND pord[1]+1 AND SO ON

  # all.equal(t(T)%*%T,diag(c1+c2-1)) # THIS IS TRUE
  # all.equal(T%*%t(T),diag(c1+c2+1)) # BUT THIS IS NOT

  Z = B%*%T_s
  X = B%*%T_n

  X_wrong=X
  T_wrong=cbind(T_n,T_s)

  X=X[,-c(pord[1],pord[1]+pord[2])] # THIS ONLY WORKS FOR pord=2 IN GENERAL IS THE SECOND COLUMN AND pord[1]+1 AND SO ON



  ####


  G=list(c(d1,rep(0,length(d2))),c(rep(0,length(d1)),d2))


  # G=list(t_1,t_2)
  # names(G)=c("t_1","t_2")
  #

  #T=rbind(cbind(T_n1,T_s1,matrix(0,nrow = dim(T_n1)[1],ncol=dim(T_n2)[2]+dim(T_s2)[2])),cbind(matrix(0,nrow = dim(T_n2)[1],ncol=dim(T_n1)[2]+dim(T_s1)[2]),T_n2,T_s2))


  ####

  list(X = X, Z = Z, G=G, T=T, T_wrong=T_wrong, X_wrong=X_wrong)
}
