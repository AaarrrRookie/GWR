
# Definition der Funktion für ein morbiSet
# Funktion entspricht dem Morbi-RSA
# angewendet auf das morbiSet nach Kost = Xb + Fehler

MRSA<-function(morbiSet,Standfehler=TRUE,Hochkosten="kein",neg_coef_control=FALSE,...){

  # morbiSet - Daten der classe morbiSet
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
  fit<-setOld(fit,morbiSet,"MRSA")       # Fülle das Objekt fit mit Informationen aus dem morbiSet

  w<-morbiSet@tage/max(morbiSet@tage)    # Definiere Gewichtung als Tage/max(Tage)

  M<-morbiSet@X@i+1                      # Definiere Ort von nicht Nullelementen von X (beachten X zählt von 0, R zählt aber von 1)


# Auswahl des Hochkostenfalls
  if(Hochkosten=="kein"){
    Y<-morbiSet@Kost/w                   # Anualisieren und Ergebnis mit dem Zeiger Y versehen
    Z_HK<-rep(0,length(Y))               # Hochkosten sind 0 und jeder hat diese Hochkosten von 0 -> Vektor nur mit 0
  }else if(Hochkosten=="KP100"){
    print("Risikopool KP100")
    HK<-KP100(morbiSet,...)              # Hochkosten aus KP100 angewendet auf morbiSet unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z                        # Hochkosten sind das Ergebnise und Ergeben extra Zuschlage
  }else if(Hochkosten=="RP1"){
    print("Risikopool RP1")
    HK<-RP1(morbiSet,...)                # Hochkosten aus KP100 angewendet auf morbiSet unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z
  }else if(Hochkosten=="RP2"){
    print("Risikopool RP2")
    HK<-RP2(morbiSet,...)                # Hochkosten aus KP100 angewendet auf morbiSet unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z
  }else if(Hochkosten=="HCP"){
    print("Risikopool HCP")
    HK<-HCP(morbiSet,...)                # Hochkosten aus KP100 angewendet auf morbiSet unter den Einstellungen ...
    Y<-HK$LA_NK/w                        # Anualisieren und Ergebnis aus Normalkosten mit dem Zeiger Y versehen
    Z_HK<-HK$HK_Z
  }else{
    stop("Keine Risikopool definiert")
  }



  print("0-Modell")

  B<-solver(morbiSet@X,Y,w,M)                         # Löse Gleichungssystem über solver

  risk.in<-c(1:ncol(morbiSet@X))                      # Alle B sind gültige Risikofaktoren (>0 und Hierarchi)

  if(neg_coef_control==TRUE){                         # wenn Verletzung überprüfen

    repeat{
      if(min(B)>0){break}                             # keine Verletzung dann stopp (break)

      print(paste0("Finde ",length(which(B<0))," kleiner Null"))

      risk.in<-risk.in[-which(B<0)]                   # Neue Risikofaktoren sind Risikofaktoren ohne Verletzung

      B<-solver(morbiSet@X[,risk.in],Y,w,M)           # Lösung durch solve
    }

    Z<-((morbiSet@X[,risk.in]%*%B)*w)@x+Z_HK          # Zuweisung sind Morbi-Zuweisung (nur geprüfte Faktoren) + Hochkosten (geg. 0)
  }else{

    Z<-((morbiSet@X%*%B)*w)@x+Z_HK                    # Zuweisung sind Morbi-Zuweisung (alle geprüfte Faktoren) + Hochkosten (geg. 0)

  }



  fit<-setZuweisung(fit,Z)                          # Erweitere das Ojekt fit um Zuweisungen


  # Prüfe die Summentreue
  # Wenn Zuweisung und Kosten mehr als 100 Euro in Summe auseinander dann Fehler

  if(abs(sum(Z-morbiSet@Kost))>100){
    warning(paste0("Fehler in der Regresion: ",sum(Z-morbiSet@Kost)))  # Sage mir das etwas Falsch ist
  }


  # Berechne Standardfehler von Wahr

if(Standfehler==TRUE){
  print("Varianz-Kovarianz")

  freiheitsG<-nrow(morbiSet@X)-(ncol(morbiSet@X))     # Freiheitsgrade sind Zeilen von X minus Spalten von X

  u<-morbiSet@Kost-Z                                  # Fehler u ist Kosten minus Zuweisung
  sigma<-(t(u)%*%u)/freiheitsG	                      # Standardabweichung der Fehler (sigma) ist Summe u^2 / Freiheit
  sigma<-sigma[1,1]                                   # sigma zu skalar (bzw. erste Element der Einelementigen Matrix)


  V<-diag(qr.solve(crossprod(morbiSet@X*w,morbiSet@X)))*sigma # Varianz der Koeffizienten siehe Standardwerke der Statistik

  fit<-setCoef(fit,V[risk.in],B)                      # Ergänze Objekt fit um Koeffizienten und deren Standardfehler

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


  Result<-morbiSet.frame(names(Z),Ausgleich_2013$Bezeichnung[M],Regression=Z,BVA=Ausgleich_2013[M,3])
  Result$Dif_abs<-Result$Regression-Result$BVA
  Result$Dif_rel<-Result$Regression/Result$BVA-1

  print(Result)

  print(sum(abs(Result$Dif_abs),na.rm=TRUE)/sum(Result$BVA,na.rm=TRUE))
  write.table(Result,paste0("Auswertung\\Regression_",stichprobe,".csv"),row.names=FALSE,sep=";",dec=",")
}

