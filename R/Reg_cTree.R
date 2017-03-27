MRSAcTree <- function(morbiSet,
                      min_Blatt=1000){


  X <- morbiSet@X
  Kost <- morbiSet@Kost
  tage <- morbiSet@tage


  Result <- data.frame(Knoten=0, Vater=NA,Name="Stamm", Kosten=sum(Kost)/sum(tage), Tage=sum(tage),Offen=0)
  Result$Name <- as.character(Result$Name)
  Edge <- rep(0,length(Kost))

  repeat{
    if(min(Result$Offen)==1){break}
    Knoten <- min(Result$Knoten[Result$Offen==0]) # Wähle den Knoten
    print(Knoten)
  i <- which(Result$Knoten==Knoten)


  Kohorte <- Edge == Knoten

  n_1 <- (tage[Kohorte]%*%X[Kohorte,])@x
  n_0 <- (Result$Tage[i]-n_1)

  m_1 <- (Kost[Kohorte]%*%X[Kohorte,])@x
  m_0 <- (Result$Kosten[i]*Result$Tage[i]-m_1)
  m_1 <- (m_1/n_1)
  m_0 <- (m_0/n_0)


  E_X <- (rep(1,length(Kost[Kohorte]))%*%X[Kohorte,]/length(X[Kohorte,]))@x
  V_K <- var(Kost[Kohorte])
  V_X <- diag(crossprod(X[Kohorte,]))-E_X^2

  COV <- Kost[Kohorte]%*%X[Kohorte,] - m_1*E_X
  COR <- COV/(sqrt(V_K)*sqrt(V_X))

  new <- which(COR==max(COR,na.rm=TRUE))

  if(length(new)==0){
    Result$Offen[i] <- 1
  }else{

    K_1 <- max(Result$Knoten)+1
    K_0 <- max(Result$Knoten)+2

    Result<-rbind(Result,data.frame(Knoten=K_1,Vater=Knoten,Name=paste0(X@Dimnames[[2]][new], " = 1"),Kosten=m_1[new],Tage=n_1[new],Offen=0))
    Result<-rbind(Result,data.frame(Knoten=K_0,Vater=Knoten,Name=paste0(X@Dimnames[[2]][new], " = 0"),Kosten=m_0[new],Tage=n_0[new],Offen=0))

    Result$Offen[i] <- 1
    Edge[X[Kohorte,new]==1] <- K_1
    Edge[X[Kohorte,new]==0] <- K_0


    Result$Offen[Result$Tage<min_Blatt] <- 1

    }

  }

}
