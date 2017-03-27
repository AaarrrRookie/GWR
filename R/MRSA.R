
# Definition der Funktion für ein theObject
# Funktion entspricht dem Morbi-RSA
# angewendet auf das theObject nach Kost = Xb + Fehler

MRSA<-function(theObject,
               Standfehler=TRUE,
               Hochkosten="kein",
               Hierarchie=FALSE,
               ...){

  # theObject - Daten der classe theObject
  # Standfehler - Standardfehler wird berechnet
  # Hochkosten - ein Hochkostenmodell wird angewendet (kein, RP1, RP2, HCP, KP100)
  # neg_coef_control - negative Regressionskoeffizienten werden eliminiert
  # ...              - weitere Einstellungen für Hochkostenmodelle (siehe jeweiliges Hochkostenmodell)


### Hilfsfunktion zur Lösung eines linearen Glleichungssystems mittels qr-Zerlegung dünnbesetzer Matrizen

  solver<-function(X,Y,w,M){

    # X - Matrix (dünn)
    # Y - Zielvektor
    # w - Gewichtungsmatrix
    # M - Matrix die den "Ort" der nicht Null Elemente von X angibt


    Xw<-X        # Zeiger Xw auf die Daten von X
    Xw@x<-w[M]   # Zeiger der Werte von Xw die nicht Null sind auf die Werte der Gewichtungsmatrix die am Gleichen Ort stehen
                 # dieser Schritt Ordnet allen Werten von Xw die nicht Null sind einen neuen Zeiger zu der auf die Werte von w zeigt,
                 # die sich ergeben würden bei X*w

    return(qr.solve(crossprod(Xw,X),crossprod(Xw,Y)))  # Lösung des Gleichungssystems (X'WX)b=X'Y mittels qr Zerlegung
  }

  fit<-fitMorbiRSA()                     # Initialisiere ein leeres Objekt vom Typ fitMorbiRSA
  fit<-setOld(fit,theObject,"MRSA")       # Fülle das Objekt fit mit Informationen aus dem theObject

  w<-theObject@tage/max(theObject@tage)    # Definiere Gewichtung als Tage/max(Tage)

  M<-theObject@X@i+1                      # Definiere Ort von nicht Nullelementen von X (beachten X zählt von 0, R zählt aber von 1)


# Auswahl des Hochkostenfalls
  if(Hochkosten=="kein"){
    Y<-theObject@Kost/w                   # Anualisieren und Ergebnis mit dem Zeiger Y versehen
    Z_HK<-rep(0,length(Y))               # Hochkosten sind 0 und jeder hat diese Hochkosten von 0 -> Vektor nur mit 0
  }else if(Hochkosten=="KP100"){
    print("Risikopool KP100")
    HK<-KP100(theObject,...)              # Hochkosten aus KP100 angewendet auf theObject unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z                        # Hochkosten sind das Ergebnise und Ergeben extra Zuschlage
  }else if(Hochkosten=="RP1"){
    print("Risikopool RP1")
    HK<-RP1(theObject,...)                # Hochkosten aus KP100 angewendet auf theObject unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z
  }else if(Hochkosten=="RP2"){
    print("Risikopool RP2")
    HK<-RP2(theObject,...)                # Hochkosten aus KP100 angewendet auf theObject unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z
  }else if(Hochkosten=="HCP"){
    print("Risikopool HCP")
    HK<-HCP(theObject,...)                # Hochkosten aus KP100 angewendet auf theObject unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z
  }else{
    stop("Keine Risikopool definiert")
  }



  print("0-Modell")

  B<-solver(theObject@X,Y,w,M)                         # Löse Gleichungssystem über solver

  if(Hierarchie==TRUE){

    data(Hierarchie_Daten)
    if(!theObject@jahr %in% Hierarchie_Daten$Jahr){
      stop(paste0("Keine Hierarchisierung für Datenjahr ",theObject@jahr," möglich"))
    }

     B_orig <- B0 <- B
     Fehler <- c()
     theObject_temp <- theObject
     repeat{
       if(sum(B<0)<1){break}
       Fehler <- c(Fehler,Null_Kontrolle(B))
       theObject_temp <- Null_Anwenden(Fehler,theObject)
       B0 <- solver(theObject_temp@X,Y,w,theObject_temp@X@i+1)

       B <- B_orig
       B[Fehler] <- 0
       B[-Fehler] <- B0

     }

     Fehler_0 <- Fehler

     B_orig <- B

     Fehler<-Hierarchie_Kontrolle(B,theObject_temp)
     theObject_temp<-Hierarchie_Anwendung(Fehler,theObject_temp)

     B0 <- solver(theObject_temp@X,Y,w,theObject_temp@X@i+1)

     B <- Hierarchie_Koef(B0,B_orig,Fehler)

     B[Fehler_0] <- 0

  }

  Z<-(theObject@X%*%B*w)@x

  fit<-setZuweisung(fit,Z)                          # Erweitere das Ojekt fit um Zuweisungen


  # Prüfe die Summentreue
  # Wenn Zuweisung und Kosten mehr als 100 Euro in Summe auseinander dann Fehler

  if(abs(sum(Z-theObject@Kost))>100){
    warning(paste0("Fehler in der Regresion: ",sum(Z-theObject@Kost)))  # Sage mir das etwas Falsch ist
  }


  # Berechne Standardfehler von Wahr

if(Standfehler==TRUE){
  print("Varianz-Kovarianz")

  freiheitsG<-nrow(theObject@X)-(ncol(theObject@X))     # Freiheitsgrade sind Zeilen von X minus Spalten von X

  u<-theObject@Kost-Z                                  # Fehler u ist Kosten minus Zuweisung
  sigma<-(t(u)%*%u)/freiheitsG	                      # Standardabweichung der Fehler (sigma) ist Summe u^2 / Freiheit
  sigma<-sigma[1,1]                                   # sigma zu skalar (bzw. erste Element der Einelementigen Matrix)


  V<-diag(qr.solve(crossprod(theObject@X*w,theObject@X)))*sigma # Varianz der Koeffizienten siehe Standardwerke der Statistik

  fit<-setCoef(fit,V,B)                      # Ergänze Objekt fit um Koeffizienten und deren Standardfehler

}
  else{
    fit<-setCoef(fit,rep(0,length(B)),B)              # Ergänze Objekt fit aber ohne Standardfehler
  }


  return(fit)                                         # Gibt fit zurück

}


## Zum Vergleich von BVA Risikofaktoren
#  aktuell nicht mehr gebraucht
Zuschlaege.MRSA<-function(fit.MRSA,GP){

  Z<-fit.MRSA$coef/365
  Z[substring(names(Z),8,10)=="AGG"]<-Z[substring(names(Z),8,10)=="AGG"]-GP/365
  nAGG<-substring(names(Z)[substring(names(Z),8,10)=="AGG"],12,14)
  NUL<-c()
  for(i in 1:length(nAGG)){NUL[i]<-paste(replicate(3-nchar(nAGG[i]),"0"),collapse="")}
  names(Z)[substring(names(Z),8,10)=="AGG"]<-paste0("AGG",NUL,nAGG)

  names(Z)[substring(names(Z),1,3)=="EMR"]<-paste0("EMG",substring(names(Z)[substring(names(Z),1,3)=="EMR"],4,7))
  names(Z)[names(Z)==""]<-"KEG001"

  Ausgleich_2013 <- read.csv2("G:/Working/1 Stichprobe/Referenzen/Ausgleich_2013.csv")

  M<-match(names(Z),Ausgleich_2013$Risikogruppe)


  Result<-theObject.frame(names(Z),Ausgleich_2013$Bezeichnung[M],Regression=Z,BVA=Ausgleich_2013[M,3])
  Result$Dif_abs<-Result$Regression-Result$BVA
  Result$Dif_rel<-Result$Regression/Result$BVA-1

  print(Result)

  print(sum(abs(Result$Dif_abs),na.rm=TRUE)/sum(Result$BVA,na.rm=TRUE))
  write.table(Result,paste0("Auswertung\\Regression_",stichprobe,".csv"),row.names=FALSE,sep=";",dec=",")
}


Null_Kontrolle <- function(B){
   for(i in 1:length(B)){
     if(B[i] < 0){
       cat(paste0(names(B)[i]," < 0 "))
       cat("\n")
     }
   }
   return(which(B<0))
}
Null_Anwenden  <- function(Fehler,theObject){
  X <- theObject@X
  X <- X[,-Fehler]
  theObject@X <- X
  return(theObject)
}

Hierarchie_Kontrolle<-function(B, theObject){
  if("HMG" %in% substring(names(B),1,3)){
    Hierarchie_Baum<-Hierarchie_Daten[Hierarchie_Daten$Jahr==theObject@jahr,c(3,4)]

    extra_null<-min(nchar(names(B)[substring(names(B),1,3)=="HMG"]))-4

    Hierarchie_Baum$A<-paste0("HMG",rep(0,extra_null),Hierarchie_Baum$A)
    Hierarchie_Baum$B<-paste0("HMG",rep(0,extra_null),Hierarchie_Baum$B)
  }else if("hmg" %in% substring(names(B),1,3)){
    Hierarchie_Baum<-Hierarchie_Daten[Hierarchie_Daten$Jahr==theObject@jahr,c(3,4)]

    extra_null<-min(nchar(names(B)[substring(names(B),1,3)=="hmg"]))-4

    Hierarchie_Baum$A<-paste0("hmg",rep(0,extra_null),Hierarchie_Baum$A)
    Hierarchie_Baum$B<-paste0("hmg",rep(0,extra_null),Hierarchie_Baum$B)
  }else{
    warnings("Keine HMG gefunden")
    return(data.frame(0))
  }

  n<-nrow(Hierarchie_Baum)
  names_B<-names(B)

  Fehler<-c()

  for(i in 1:nrow(Hierarchie_Baum)){
      B_A<-B[which(names_B==Hierarchie_Baum$A[i])]
      B_B<-B[which(names_B==Hierarchie_Baum$B[i])]

      if(length(B_A)>0&length(B_B)>0){
        if(B_B>B_A){
          cat("Hierarchieverletzung: ")
          cat(paste0(Hierarchie_Baum$A[i]," < ",Hierarchie_Baum$B[i]))
          cat("\n")
          Fehler<-c(Fehler,i)
        }
      }
  }

  if(length(Fehler)==0){
    return(data.frame(0))
  }else{
    Fehler<-Hierarchie_Baum[Fehler,]
    return(Fehler)
  }

}

Hierarchie_Anwendung<-function(Fehler,theObject){
  if(nrow(Fehler)==1){
    return(theObject)
  }else{
    X<-theObject@X
    names_X<-X@Dimnames[[2]]

    for(i in 1:nrow(Fehler)){
      x_A<-which(names_X==Fehler$A[i])
      x_B<-which(names_X==Fehler$B[i])

      if(length(x_A)>0&length(x_B)>0){

        X[,x_A]<-X[,x_A]+X[,x_B]
        X <- X[,-x_B]

        names_X<-X@Dimnames[[2]]
      }
    }
  }

  theObject@X<-X

  return(theObject)

}


Hierarchie_Koef <- function(B0,B_orig,Fehler){

  B <- B0
  n_B0 <- names(B0)
  n_B  <- names(B_orig)

  for(i in 1:length(B_orig)){
    if(n_B[i] %in% n_B0){
       a <- n_B[i]
       B[i] <- B0[which(n_B0==n_B[i])]
    }else{
       a <- n_B[i]
       aa <- Fehler$A[Fehler$B==a][1]
       if(aa %in% n_B0){
         B[i] <- B0[which(n_B0==aa)]
       }
    }
  }

  names(B) <- names(B_orig)

  return(B)
}
