#### Funktion für die Berechnung des Morans'I nach Durchführung einer Regression (beliebig)



MoranI.Test<-function(Data,fit,k_krit=5000,Kernel="bi",typ="normal",n_sim=100){
  # Data - Daten vom Typ morbiSet
  # fit  - Ergebnis einer Regression
  # k_krit - Grenzwert für die Distanz / Versichertenanzahl (Standard 5000)
  #!!! Bitte Anpassen von k_krit auf Kernel und sinnvolle Werte (mehr als einen?) -> Sonst sinnlose Ergebnisse
  #Kernel - Distanzgewichtungsfunktion (bi, inf, GWR)
  #typ - Annahme an die Verteilung des MI (normal, MC) als normalverteilt oder unbekannt (Monte Carlo)
  #n_sim - Anzahl Simulationen

  data(OSM)

  IN<-Data@plz %in% OSM$plz				  # Nur Versicherte mit bekanntem Wohnort verwenden

  plz<-as.integer(Data@plz[IN])     # Postleitzahlen zu integer

  print("Distanz")

  Dist<-distanz(OSM,Data@plz[IN])   # Berechne Distanzmatix

  print("Gewichtung")


  # Auswahl und Berechnung der Gewichtungsfunktion
  if(Kernel=="bi"){
    W<-WEIGTH.bi(Dist$D,k_krit)					# Aus Distanz folgende Gewichtungsamtrix
  }else if(Kernel=="inf"){
    W<-WEIGTH.inf(Dist$D,k_krit)
  }else if(Kernel=="GWR"){
    W<-WEIGTH(Dist,k_krit)
  }else{
    stop("Falscher Kern")
  }


  E<-(fit@Zuweisung-fit@Kost)         # Regressionsfehler
  TSS<-sum(E[IN]^2)                   # total sum of squares / zweites Moment der Verteilung der Fehler
  y4<-sum(E[IN]^4)                    # viertes Moment der Verteilung der Fehler

  I<-MoranI(as.integer(Dist$Match-1),plz[IN],Dist$Dplz,W,E[IN],TSS)    # Morans'I Funktion siehe Funktionen in "C"
  n<-I$n                              # Anzahl Beobachtung für MI

  EI<-(-1/(I$n-1))                    # Erwartungswert von MI
  S1<-I$s1                            # S1 (Stützwert für Varianz)
  S2<-4*n                             # S2 (Stützwert für Varianz)
  S3<-(y4/n)/(TSS/n)^2                # S3 (Stützwert für Varianz)
  S4<-(n^2-3*n+3)*S1-n*S2+3*n^2       # S4 (Stützwert für Varianz)
  S5<-(n^2-n)*S1-2*n*S2+6*n^2         # S5 (Stützwert für Varianz)

  V<-(n*S4-S3*S5)/((n-1)*(n-2)*(n-3)*n^2)-EI^2  # Varianz des MI unter Normalverteilungsannahme

  if(typ=="MC"){
    print("Monte Carlo Simulation")
    print("!ACHTUNG, sehr lange Laufzeit")
    I_MC<-rep(NA,n_sim)
    for(i in 1:n_sim){
      I_MC[i]<-MoranI_MC(sample(as.integer(Dist$Match-1)),plz[IN],Dist$Dplz,W,E[IN],TSS)    # Morans'I Funktion siehe Funktionen in "C"
      if(i%%n_sim==0){
        print(i/n_sim)
      }
    }
    EI<-mean(I_MC)                                             # Erwartungswert
    V<-var(I_MC)                                               # Varianz
  }


  MIZ<-(I$I-EI)/sqrt(V)                                       # MI Standardisiert
  p<-(1-pt(MIZ,length(Data@tage)))*2                          # p-Wert

  return(list(I=I,MIZ=MIZ,p=p))                               # Rückgabe der Werte


}



#Aktuell nicht implementiert
# Eingesetzt in der Publikation Gesundheitswesen aktuell
test.region<-function(Data,fit){

  Ergebnis<-list()

  plz<-as.integer(as.character(Data$plz))
  OSM <- read.csv("G:/Regionaldaten/OSM_P.csv")
  IN<-plz %in% OSM$plz

  Dist<-distanz(OSM,plz[IN])

  print("MI50")
  #### MI50

  W<-WEIGTH.bi(Dist$D,50000)
  TSS<-sum(fit$E[IN]^2)
  y4<-sum(fit$E[IN]^4)

  I<-MoranI(as.integer(Dist$Match-1),plz[IN],Dist$Dplz,W,fit$E[IN],TSS)
  n<-I$n

  EI<-(-1/(I$n-1))
  S1<-I$s1
  S2<-4*n
  S3<-(y4/n)/(TSS/n)^2
  S4<-(n^2-3*n+3)*S1-n*S2+3*n^2
  S5<-(n^2-n)*S1-2*n*S2+6*n^2

  V<-(n*S4-S3*S5)/((n-1)*(n-2)*(n-3)*n^2)-EI^2

  MIZ<-(I$I-EI)/sqrt(V)
  p<-(1-pt(MIZ,nrow(Data)))*2


  Ergebnis$I.50=list(I=I$I,MIZ=MIZ,p=p)


  ### MI80
  print("MI80")

  W<-WEIGTH.bi(Dist$D,80000)
  TSS<-sum(fit$E[IN]^2)
  y4<-sum(fit$E[IN]^4)

  I<-MoranI(as.integer(Dist$Match-1),plz[IN],Dist$Dplz,W,fit$E[IN],TSS)
  n<-I$n

  EI<-(-1/(I$n-1))
  S1<-I$s1
  S2<-4*n
  S3<-(y4/n)/(TSS/n)^2
  S4<-(n^2-3*n+3)*S1-n*S2+3*n^2
  S5<-(n^2-n)*S1-2*n*S2+6*n^2

  V<-(n*S4-S3*S5)/((n-1)*(n-2)*(n-3)*n^2)-EI^2

  MIZ<-(I$I-EI)/sqrt(V)
  p<-(1-pt(MIZ,nrow(Data)))*2

  Ergebnis$I.80=list(I=I$I,MIZ=MIZ,p=p)

  print("LME")
  ### LME

  LME.value<-LME(as.integer(Dist$Match-1),plz[IN],Dist$Dplz,W,fit$E[IN],TSS)

  p<-1-pchisq(LME.value$T,1)


  Ergebnis$LME=list(T=LME.value$T,p)

  print("LML")
  ### LML

  sigma<-LME.value$s
  T<-LME.value$s1

  BWY<-GWR_core1(as.integer(Dist$Match-1),W,fit$Kost[IN])
  B.WYhat<-GWR_core1(Dist$Match,W,fit$Z[IN])
  WYhat<-B.WYhat[Dist$Match]
  WY<-BWY[Dist$Match]

  WSS<-sum(fit$E[IN]*WY)

  fit.WYHAT<-qr.solve(crossprod(X[IN,],X[IN,]),crossprod(X[IN,],WYhat))

  T2<-sum(WYhat*(WYhat-(X[IN,]%*%fit.WYHAT)@x))

  T<-T+T2/TSS

  LML.value<-(WSS/TSS)^2/T

  p<-1-pchisq(LML.value,1)


  Ergebnis$LML=list(T=LML.value,p)

  return(Ergebnis)

}

# Funktionen für einen RESET-Test
# Implementierung später geplant
test.reset<-function(fit){

  library(speedglm)

  Y<-fit$Kost
  Yhat<-fit$Z

  TSS<-sum((Y-mean(Y))^2)

  X<-matrix(Yhat/mean(Yhat))
  #X<-cbind(X,1)

  R0<-1-sum((Y-X%*%qr.solve(crossprod(X,X),crossprod(X,Y)))^2)/TSS

  X<-cbind(X,X[,1]^2)
  X<-cbind(X,X[,1]^3)

  R1<-1-sum((Y-X%*%qr.solve(crossprod(X,X),crossprod(X,Y)))^2)/TSS



  n<-length(Y)
  m<-length(fit$coef)
  k<-3

  T<-((R1-R0)/k)/((1-R1)/(n-(n-m)-1))
  p<-1-pf(T,k-1,n-m-k)

}
