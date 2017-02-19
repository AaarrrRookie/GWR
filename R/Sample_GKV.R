
Sample <- function(X,Ref,n){

  R <- Ref$Anteil[match(X@Dimnames[[2]],Ref$Risikogruppe)]
  R <- R*n

  P <- 1/2-1/2*X%*%solve(crossprod(X))%*%(t(X)%*%rep(1/2,nrow(X))-2*R)

  return(P)

}
