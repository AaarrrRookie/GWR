distanz<-function(plz_P,plz){

  Coord<-plz_P
  Dplz<-Coord$plz
  D<-matrix(NA,length(Dplz),length(Dplz))
  for(i in 1:length(Dplz)){
    x<-rep(Coord[i,2],length(Dplz))-Coord[,2]
    y<-rep(Coord[i,3],length(Dplz))-Coord[,3]
    D[,i]<-sqrt(x^2+y^2)
  }
  remove(x,y)
  Match<-match(plz,Dplz)

  return(list(D=D,Coord=Coord,Match=Match,plz=plz,Dplz=Dplz))
}


distanz_zeit<-function(plz_P,plz,zeit){

  Coord<-plz_P
  Dplz<-Coord$plz
  D_O<-matrix(NA,length(Dplz),length(Dplz))
  for(i in 1:length(Dplz)){
    x<-rep(Coord[i,2],length(Dplz))-Coord[,2]
    y<-rep(Coord[i,3],length(Dplz))-Coord[,3]
    D_O[,i]<-sqrt(x^2+y^2)
  }


  D_zeit <- unique(zeit)
  D_T <- matrix(NA,length(D_zeit),length(D_zeit))

  for(i in 1:length(D_zeit)){
    x<-rep(D_zeit[i],length(D_zeit))-D_zeit
    D_T[,i] <- abs(x)
  }


  remove(x,y)
  Match<-match(plz,Dplz)
  Match_T <- match(zeit,D_zeit)

  return(list(D=D_O,Match=Match,plz=plz,Dplz=Dplz,
              D_T = D_T, Match_T=Match_T , zeit=zeit, D_zeit=D_zeit))
}




WEIGTH<-function(D,k){
  W<-matrix(0,nrow(D$D),ncol(D$D))
  n<-length(D$Match)

  Vers<-data.table(Vers=rep(1,n),plz=D$plz)
  Vers<-data.frame(Vers[,sum(Vers),by=plz])
  M<-match(D$Dplz,Vers$plz)
  V<-Vers[M,]
  V$V1[is.na(V$V1)]<-0

  for(i in 1:ncol(D$D)){
    o<-order(D$D[,i])
    A<-V[o,]
    K<-which(cumsum(A$V1)<k)
    if(length(K)==0){
      K<-1
    }
    lK<-K[length(K)]
    g<-D$D[o[lK],i]
    W[o[1:lK],i]<-(1-D$D[o[1:lK],i]^2/g^2)^2
  }
  diag(W)<-1
  return(W)
}





WEIGTH.zeit<-function(D,k){
  W<-matrix(0,nrow(D$D_T),ncol(D$D_T))
  n<-length(D$Match_T)

  Vers<-data.table(Vers=rep(1,n),zeit=D$zeit)
  Vers<-data.frame(Vers[,sum(Vers),by=zeit])
  M<-match(D$D_zeit,Vers$zeit)
  V<-Vers[M,]
  V$V1[is.na(V$V1)]<-0

  for(i in 1:ncol(D$D_T)){
    o<-order(D$D[,i])
    A<-V[o,]
    K<-which(cumsum(A$V1)<k)
    if(length(K)==0){
      K<-1
    }
    lK<-K[length(K)]
    g<-D$D[o[lK],i]
    W[o[1:lK],i]<-(1-D$D[o[1:lK],i]^2/g^2)^2
  }
  diag(W)<-1
  return(W)
}




WEIGTH.bi<-function(D,k){
  W<-(1-(D/k)^2)^2
  W[D>k]<-0

  return(W)
}

WEIGTH.inv<-function(D,k){
  W<-1/D^k

  return(W)
}
