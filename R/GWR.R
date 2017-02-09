# Definition der Funktion für ein morbiSet
# Funktion entspricht dem Morbi-RSA mit GWR-Regionalfaktor
# angewendet auf das morbiSet nach Kost = Xb + GWR + Fehler

GWR<-function(morbiSet,k_Krit,exakt=FALSE,std_exakt=FALSE,treu=FALSE,Kernel="GWR",...){	# GWR-Regression

  # morbiSet - Daten der classe morbiSet
  # k_krit   - kritischer Wert für die Gewichtung
  # exakt    - GWR exakt oder aproximativ (exakt nur bei kleinem Datensatz notwenig)
  # Standfehler - Standardfehler wird berechnet
  # Hochkosten - ein Hochkostenmodell wird angewendet (kein, RP1, RP2, HCP, KP100)
  # neg_coef_control - negative Regressionskoeffizienten werden eliminiert
  # ...              - weitere Einstellungen für Hochkostenmodelle (siehe jeweiliges Hochkostenmodell)



  if(k_Krit<0){
    stop("Kritischer Wert nicht angegeben")
  }

 data(OSM)                           # Lade Georeferenzen aus OSM (liegt der Bibliothek bei)

  print("Auslandsversicherte")

  IN<-morbiSet@plz %in% OSM$plz				# Nur Versicherte mit bekanntem Wohnort verwenden, alle Anderen als Wahr in dem Vektor IN markieren

  plz<-as.integer(morbiSet@plz[IN])   # Kopiere in plz von denjenigen die bekannt sind -< als integer

  print("Distanz")

  Dist<-distanz(OSM,morbiSet@plz[IN])			          # Distanz aller Versicherten zueinander mittels distanzfunktion

  print("Gewichtung")

  # Auswahl und Berechnung der Gewichtungsfunktion
  if(Kernel=="bi"){
    W<-WEIGTH.bi(Dist$D,k_Krit)					# Aus Distanz folgende Gewichtungsamtrix
  }else if(Kernel=="inf"){
    W<-WEIGTH.inv(Dist$D,k_Krit)
  }else if(Kernel=="GWR"){
    W<-WEIGTH(Dist,k_Krit)
  }else{
    stop("Falscher Kern")
  }

  print("MRSA")

  if(sum(IN)<length(morbiSet@Kost)){
    warning("Versicherte ohne PLZ gefunden")
    morbiSet@X<-cBind(morbiSet@X,!IN)  # Wenn es unbekannte PLz gibt, dann sind diese als Restgröße (Dummy) an die Designmatrix anzuhängen (cBind kopiert links dran)
    M<-match(morbiSet@plz,Dist$Dplz)
  }else{
    M<-match(morbiSet@plz,Dist$Dplz)									# Rückzuordnung der GWR-Koeffizienten zu Versicherten über PLZ
  }
  if(exakt==FALSE){
    fit.MRSA<-MRSA(morbiSet,...)	     # nicht räumliches MRSA mit Auslandsfaktor (ungekannte PLZ)
    E<-fit.MRSA@Kost-fit.MRSA@Zuweisung							# E - Deckungslücke des MRSA
  }



  ### Exakte GWR-Regression (nur bei kleinem Datensatz)
  # Anpassung des X auf Korrelation zur Regionsvariable
  if(exakt==TRUE){

    #X<-as.matrix(morbiSet@X)
    #storage.mode(X)<-"integer"                                # Wandle X zurück in dicke Matrix

    X<-t(morbiSet@X)

    Sv<-GWR_core2(as.integer(Dist$Match-1),W,X@i,X@p,length(morbiSet@Kost),ncol(morbiSet@X)) # GWR für X zur Berechnung der räumlichen Effekte in X

    #X<-morbiSet@X-Sv[M,]                                     # X bereinigt um räumliche Effekte von X



    w<-morbiSet@tage/max(morbiSet@tage)                       # Gewichtungsvektor für Versichtentage
    Y<-morbiSet@Kost/w                                        # Anualisierung (Kosten durch Gewichtung)

    XX<-as((X*w)%*%morbiSet@X,"matrix")

    SS<-cross_vec(as.integer(Dist$Match-1),Sv,length(morbiSet@Kost),ncol(morbiSet@X))
    SX<-cross_mat(as.integer(Dist$Match-1),X@i,X@p,Sv,length(morbiSet@Kost),ncol(morbiSet@X))

    VarCov<-XX-2*SX+SS
    Z0<-cross_mat_vec(as.integer(Dist$Match-1),Sv,Y,length(morbiSet@Kost),ncol(morbiSet@X))
    Z<-((X*w)%*%Y-Z0)@x

    B<-qr.solve(VarCov,Z)                                     # Lösen des linearen Gleichungssystems (x'WX)b=X'WY
    E0<-cross_vec_mat(as.integer(Dist$Match-1),Sv,B,length(morbiSet@Kost),ncol(morbiSet@X))
    E<-(morbiSet@Kost-morbiSet@X%*%B*w)@x-E0      # Fehler der Regression
  }



#################################################################################################
#                                                                                               #
#                                                                                               #
#                           GWR-Regression				                                        #
#                                                                                               #
#################################################################################################
    print("GWR")

  if(exakt==FALSE){
    coef.REG<-GWR_core1(as.integer(Dist$Match-1),W,E[IN]-mean(E[IN]))		                      # Regressionskoeffizienten aus GWR
    names(coef.REG)<-Dist$Dplz                                                        # Benenne Koeffizienten nach PLZ
  }else{
    coef.REG<-GWR_core1(as.integer(Dist$Match-1),W,E[IN])   		                      # Regressionskoeffizienten aus GWR
    names(coef.REG)<-Dist$Dplz                                                        # Benenne Koeffizienten nach PLZ
  }
   # GWR<-data.frame(plz=Dist$Dplz,Gruppe=Dist$Dplz,COEF=coef.REG)
   # GWR$Gruppe<-as.character(GWR$Gruppe)												# PLZ wieder zu String
   # GWR$Gruppe[nchar(GWR$Gruppe)==4]<-paste0("0",GWR$Gruppe[nchar(GWR$Gruppe)==4])


   # morbiSet.REG<-data.frame(Kosten=morbiSet@Kost,Z=fit.MRSA@Zuweisung,REG=coef.REG[M])  # definiere Data.frame aus morbiSet und Ergebnissinformationen
   #  morbiSet.REG$REG[is.na(morbiSet.REG$REG)]<-0    # Dataframe ohne Regionalzuschlag hat den Zuschlag 0

    Z.Reg<-(coef.REG[M]/365)*morbiSet@tage
    Z.Reg[is.na(Z.Reg)]<-0

    GP.r<-mean(Z.Reg)           										# Anpassung der Grundpauschale an GWR-Koeffizienten (Null-Effekt über Region)
    coef.REG<-coef.REG-GP.r                         # Anpassen der Koeffizienten
    Z.Reg<-(coef.REG[M]/365)*morbiSet@tage
    Z.Reg[is.na(Z.Reg)]<-0

   # morbiSet.REG$E<-morbiSet.REG$Kosten-morbiSet.REG$Z-morbiSet.REG$REG # Fehler der Regression (Kosten - fixe Zuweisung - GWR)

    if(treu==TRUE){
      morbiSet@Kost<-morbiSet@Kost-Z.Reg
      fit.MRSA<-MRSA(morbiSet,...)	     # nicht räumliches MRSA mit Auslandsfaktor (ungekannte PLZ)
      fit.MRSA@Kost<-fit.MRSA@Kost+Z.Reg
    }

  Z<-fit.MRSA@Zuweisung+Z.Reg               # Zuweisung fix + GWR
  E<-morbiSet@Kost-Z                        # Fehler

  fit<-setZuweisung(fit.MRSA,Z)                  # Fülle Ergebnisobjekt mit Informationen


#################################################################################################
#                                                                                               #
#                                                                                               #
#                          Varianz-Kovarianz			                                        #
#                                                                                               #
#################################################################################################



  ### COEF
  print("Varianz - Kovarianz")									#


  sigma<-sd(E)                                  # Standardabwichung der Fehler

  if(exakt==FALSE){
    B.MRSA<-fit.MRSA@Coef$B                     # wenn exakt dann B sonst B aus Ergebnisobjekt vom Morbi-RSA
    #B.MRSA[1:40]<-B.MRSA[1:40]-GP.r             # Anpassen der Grundpauschale (hier AGG 1 bis 40) an Regionaleffekt
    names(B.MRSA)<-rownames(fit.MRSA@Coef)      # jeweils richtige Namen zuordnen
  }else{
    B.MRSA<-B                                   # ...
    #B.MRSA[1:40]<-B.MRSA[1:40]-GP.r
    names(B.MRSA)<-names(B)                     # ...
  }

 # SD.MRSA<-diag(qr.solve(crossprod(morbiSet@X*morbiSet@tage/365,morbiSet@X)))*sigma   # Standardabweichung fixe Koeffizienten siehe Standardlehrbuch
   SD.MRSA<-fit.MRSA@Coef$SD^2


  ### GWR
  B.GWR<-coef.REG

  if(std_exakt==TRUE){
    coef.V<-GWR_core1(as.integer(Dist$Match-1),W,E^2)
    SD.GWR<-coef.V

  }else{


    Vers<-data.table(plz=Dist$plz,Vers=1)                    # Datensatz bestehend aus PLZ und immer 1
    Vers<-Vers[,sum(Vers),by=plz]                            # Summiere über PLZ die data.table Art (viel schneller als apply)
    Vers<-Vers[match(Dist$Dplz,Vers$plz),2]                  # Ordne (match) den Regressionskoeffizienten die Anzahl an Versicherten zu
    V<-colSums(W^2*Vers$V1,na.rm=TRUE)                          # Schritt 1: Die Summe der Quadrierten GWR-Gewichte
    V<-V^2/Vers/k_Krit                                       # Schritt 2: Quadrieren durch Anzahl der Versicherten und kritischen Wert Teilen
    SD.GWR<-c(V$V1*sigma)^2                                     # V*sigma ist Standardabwichung der GWR-Koeffizienten

  }

  fit<-setCoef(fit,c(SD.MRSA,SD.GWR),c(B.MRSA,B.GWR))  # füttere Ergebnisobjekt mit Informationen


  return(fit)  # Rückgabe Ergebnisobjekt

}


# Hilfsfunktion die nur GWR-Koeffizienten ohne Test, ohne MRSA und ohne Varianz erzeugt (schnell)
GWR_link<-function(Y,k_Krit,plzVers){

  # Y 		= Zielvariable
  # k_Krit	= Erwartete Bandbreite (50 bis 100 * 1000 Versicherte)
  # plzVers = PLZ der Daten
  # mit morbiSet starten als GWR_link(morbiSet@Kost, z.B. 10000, morbiSet@plz)

  # Rest wie oben (unkommentiert)



  print("Auslandsversicherte")

  IN<-plzVers %in% OSM$plz

  plz<-as.integer(plzVers[IN])

  print("Distanz")

  Dist<-distanz(OSM,plzVers[IN])

  print("Gewichtung")

  W<-WEIGTH(Dist,k_Krit)

  print("GWR")

  coef.REG<-GWR_core1(as.integer(Dist$Match-1),W,as.numeric(Y))
  names(coef.REG)<-Dist$Dplz

  return(coef.REG)

}



# Cross-Validierung über GWR
CV_GWR<-function(morbiSet,part=10,chain=1,k=seq(50,100,10),...){		# Kreuzvalidierung über GWR

  # morbiSet = Daten als morbiSet
  # part	= Größe und Anzahl der Testmengen (10 - 10 Testmengen a 1/10 gelernt auf den jeweils anderen 9/10)
  # chain = wieviele Zufallsstarts
  # k		= Erwarteter Bereich der Bandbreite (50 bis 100 * 1000 Versicherte)

  data(OSM)   # Lade Referenzsdaten aus OSM


  IN<-morbiSet@plz %in% OSM$plz				# Nur Versicherte mit bekanntem Wohnort verwenden, alle Anderen als Wahr in dem Vektor IN markieren
  morbiSet<-subsetMorbiSet(morbiSet,IN)

  k<-k*1000

  # Definiere Prüfgrößen

  R2<-MAPE<-CPM<-matrix(NA,(chain),length(k))     # Container für Validierungsmenge hier R2, MAPE und CPM
  N<-POOL<-CV<-c(1:length(k))                     # Durchlaufparameter

  # Zufällige Testmengen

  set.seed(1)                                     # Fester Wert für Zufallsstart

  n_morbiSet<-nrow(morbiSet@X)                    # Wahre Anzahl an Beobachtungen
  frac<-rep(c(1:part),n_morbiSet/part)            # Zerlege Morbiset in Teile definierter größe
  frac<-sample(frac)                              # Ziehe zufällig aus den Teilen


  for (j in 1:chain){# für Schritt 1 bis chain


    for (i in 1:length(k)){ # für Teil 1 bis k

      print(paste("k= ",k[i]," ON CHAIN: ",j,sep=""))

      plz<-morbiSet@plz[!frac==j]

      subset<-frac==j                                     # Definiere Untermenge vom morbiSet
      morbiSet.test<-subsetMorbiSet(morbiSet,subset)     # morbiSet Testmenge
      morbiSet.lern<-subsetMorbiSet(morbiSet,!subset)      # morbiSet Lernmenge


      fit.MRSA<-MRSA(morbiSet.lern,Standfehler=FALSE,...)			# MRSA auf Lernmenge

      E<-fit.MRSA@Kost-fit.MRSA@Zuweisung                 # Fehler auf Lernmenge

      Dist<-distanz(OSM,plz)														  # Distanz der Lernmenge
      W<-WEIGTH(Dist,k[i]*(1-1/part))											# GWR-Gewicht auf Lernmenge


      print("GWR")

      coef.REG<-GWR_core1(as.integer(Dist$Match-1),W,E-mean(E))		# Regressionskoeffizienten

      ############################# Auswerten der Ergebnisse auf Testmenge

      Z.MRSA<-(morbiSet.test@X%*%fit.MRSA@Coef$B)@x/365  # Zuweisung fix auf Testmenge
      M<-match(morbiSet.test@plz,Dist$Dplz)              # Zuordnung Testmengen plz zu Ergebnis-plz der Lernmenge (M sagt welches Element zu welchem)
      Z.GWR<-coef.REG[M]/365                             # Zuordnung der GWR-Koeffizienten mittels Zuordnungsmenge M


      Z<-(Z.MRSA+Z.GWR)*morbiSet.test@tage               # Zuweisung gesamt aus fix + morbiRSA
      K<-morbiSet.test@Kost                              # Kosten der Testmenge


      # Berechne und Schreibe Indikatoren R2, Mape, CPM
      R2[j,i]<-var(Z)/var(K)
      MAPE[j,i]<-sum(abs(K-Z))/length(K)
      CPM[j,i]<-1-sum(abs(K-Z))/sum(abs(K-mean(K)))
      print(MAPE[j,i])
    }

  }

  CV<-colSums(MAPE)/ncol(MAPE)                              # Kriterium ist R2
  plot(k,CV)                                            # plotte (Grafik) die Ergebnisse


  return(list(CV=CV,R2=R2,MAPE=MAPE,CPM=CPM,k=k))       # Rückgabe der Ergebnisse als Liste
}
