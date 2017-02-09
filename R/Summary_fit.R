######################################################################################################

#                                            Definiere Standard-Prozeduren

######################################################################################################


# Prozedur für Summary
# Wenn summary(fit)
# gibt die Koeffizientenmatrix mit Standardabweichung und t-Test aus
# Sonst keine Funktionen integriert (jetzt)
summary.fitMorbiRSA<-function(theObject){
  print(paste0(theObject@Modell," mit ",length(theObject@Kost)," Beobachtungen"))
  print(theObject@Coef)
  print(theObject@Guete)
}


write.fit<-function(morbiSet
                    ,fit.MRSA=fitMorbiRSA()   # leeres Objekt vom Typ fit
                    ,fit.Reg=fitMorbiRSA()    # leeres Objekt vom Typ fit
                    ,file=""                  # Name der Ergebnisdatei (nur .xlsx)
){

  wb <- createWorkbook()                                  # generiere Arbeitsblatt
  options("openxlsx.borderStyle" = "thin")                # Tabellen-Definition
  options("openxlsx.borderColour" = "#4F81BD")

  n<-length(morbiSet@Kost)                   # Anzahl der Daten
  n.MRSA<-length(fit.MRSA@Kost)              # Anzahl der Daten MRSA
  n.Reg<-length(fit.Reg@Kost)                # Anzahl der Daten Reg


  # Auswahl des Ergebnistyp
  if(!n==n.MRSA&!n==n.Reg){
    stop("Ergebnisse passen nicht zum MorbiSet")
  }else if(n==n.MRSA&n==n.Reg){
    print("MRSA-Vergleich")
    type<-"vergleich"
  }else if(n==n.Reg){
    print("Auswertung Regionalmodell")
    type<-"nur_Regio"
  }else{
    print("Auswertung Morbi-RSA")
    type<-"nur_MRSA"
  }

  ##### Koeffizienten schreiben nach Typ

  if(type=="vergleich"){
    addWorksheet(wb, "COEF_MRSA")                               # Generiere Arbeitsblatt
    addWorksheet(wb, "COEF_Reg")

    writeData(wb, "COEF_MRSA", fit.MRSA@Coef, rowNames = TRUE)  # Fülle Arbeitsblatt
    writeData(wb, "COEF_Reg", fit.Reg@Coef, rowNames = TRUE)    # Fülle Arbeitsblatt
  }else if(type=="nur_MRSA"){
    addWorksheet(wb, "COEF_MRSA")

    writeData(wb, "COEF_MRSA", fit.MRSA@Coef, rowNames = TRUE)
  }else{
    addWorksheet(wb, "COEF_Reg")

    writeData(wb, "COEF_Reg", fit.Reg@Coef, rowNames = TRUE)
  }


  ############################## Guete schreiben nach Typ siehe Oben

  if(type=="vergleich"){
    addWorksheet(wb, "Guete")

    Guete<-fit.MRSA@Guete          # data.frame "Guete" aus MRSA
    Guete[2,]<-fit.Reg@Guete       # data.frame "Guete" aus Regionalmodell
    Guete<-data.frame(Modell=c("MRSA","Reg"),Guete)  # Zusammenführen und benennen (erste Spalte hat Namen)

    writeData(wb, "Guete", Guete, rowNames = FALSE)  # Fülle Arbeitsblatt

  }else if(type=="nur_MRSA"){
    addWorksheet(wb, "Guete_MRSA")

    writeData(wb, "Guete_MRSA", fit.MRSA@Guete, rowNames = TRUE)
  }else{
    addWorksheet(wb, "Guete_Reg")

    writeData(wb, "Guete_Reg", fit.Reg@Guete, rowNames = TRUE)
  }


  ############################### Tablen erstellen

  LA<-morbiSet@Kost               # Kosten aus dem morbiSet (Zeiger)

  # Zeiger auf Zuweisung wenn Vorhanden sonst NA
  if(type=="vergleich"){
    Z.MRSA<-fit.MRSA@Zuweisung
    Z.Reg<-fit.Reg@Zuweisung
  }else if(type=="nur_MRSA"){
    Z.MRSA<-fit.MRSA@Zuweisung
    Z.Reg<-rep(NA,length(Z.MRSA))
  }else{
    Z.Reg<-fit.Reg@Zuweisung
    Z.MRSA<-rep(NA,length(Z.Reg))
  }


  ##### AGG - Arbeitsblatt

  AGG<-morbiSet@X[,1:40]            # AGG ist eine Matrix (0/1) wenn AGG zutrifft oder nicht aus Designmatrix in morbiSet
  addWorksheet(wb, "AGG")

  Table<-data.frame(Kosten=(LA%*%AGG)@x,Zuw.MRSA=(Z.MRSA%*%AGG)@x,Zuw.Reg=(Z.Reg%*%AGG)@x,N=(rep(1,nrow(AGG))%*%AGG)@x)    # Summen der Spalten der AGG gibt absolute Summen je AGG
  Table$Kosten<-Table$Kosten/Table$N            # Teilen durch Versichertenanzahl ja AGG
  Table$Zuw.MRSA<-Table$Zuw.MRSA/Table$N        # ...
  Table$Zuw.Reg<-Table$Zuw.Reg/Table$N
  Table$Dif<-Table$Zuw.Reg-Table$Zuw.MRSA       # Differenz der Modelle
  Table$Deck.MRSA<-Table$Zuw.MRSA/Table$Kosten  # Deckungsquoten (Zuweisung / Kosten) je AGG
  Table$Deck.Reg<-Table$Zuw.Reg/Table$Kosten
  Table$Versicherte<-Table$N                    # Nenne "N" ab jetzt "Versicherte"
  Table$N<-NULL                                 # Lösche Spalte N (vorher "Versicherte" genannt)

  writeData(wb, "AGG", Table, rowNames = TRUE)  # Füllen des Blattes AGG


  ##### EMG

  IN_EMG<-which(substring(colnames(morbiSet@X),1,3)=="EMG")  # Test ob und wo der String "EMG" auftaucht

  if(length(IN_EMG)>0){       # Nur wenn EMG vorhanden
    EMG<-morbiSet@X[,IN_EMG]  # Stelle in der Designmatrix mit EMG als neue MAtrix mit Namen EMG
    # Rest wie AGG
    addWorksheet(wb, "EMG")

    Table<-data.frame(Kosten=(LA%*%EMG)@x,Zuw.MRSA=(Z.MRSA%*%EMG)@x,Zuw.Reg=(Z.Reg%*%EMG)@x,N=(rep(1,nrow(AGG))%*%EMG)@x)
    Table$Kosten<-Table$Kosten/Table$N
    Table$Zuw.MRSA<-Table$Zuw.MRSA/Table$N
    Table$Zuw.Reg<-Table$Zuw.Reg/Table$N
    Table$Dif<-Table$Zuw.Reg-Table$Zuw.MRSA
    Table$Deck.MRSA<-Table$Zuw.MRSA/Table$Kosten
    Table$Deck.Reg<-Table$Zuw.Reg/Table$Kosten
    Table$Versicherte<-Table$N
    Table$N<-NULL

    writeData(wb, "EMG", Table, rowNames = TRUE)
  }


  ##### KEG - wie Oben

  IN_KEG<-which(substring(colnames(morbiSet@X),1,3)=="KEG")

  if(length(IN_KEG)>0){
    KEG<-morbiSet@X[,IN_KEG]

    addWorksheet(wb, "KEG")

    Table<-data.frame(Kosten=(LA%*%KEG)@x,Zuw.MRSA=(Z.MRSA%*%KEG)@x,Zuw.Reg=(Z.Reg%*%KEG)@x,N=colSums(KEG))
    Table$Kosten<-Table$Kosten/Table$N
    Table$Zuw.MRSA<-Table$Zuw.MRSA/Table$N
    Table$Zuw.Reg<-Table$Zuw.Reg/Table$N
    Table$Dif<-Table$Zuw.Reg-Table$Zuw.MRSA
    Table$Deck.MRSA<-Table$Zuw.MRSA/Table$Kosten
    Table$Deck.Reg<-Table$Zuw.Reg/Table$Kosten
    Table$Versicherte<-Table$N
    Table$N<-NULL

    writeData(wb, "KEG", Table, rowNames = TRUE)

  }
  ##### HMG - wie Oben

  addWorksheet(wb, "HMG")

  HMG<-morbiSet@X[,41:ncol(morbiSet@X)]

  Table<-data.frame(Kosten=(LA%*%HMG)@x,Zuw.MRSA=(Z.MRSA%*%HMG)@x,Zuw.Reg=(Z.Reg%*%HMG)@x,N=(rep(1,nrow(AGG))%*%HMG)@x)
  Table$Kosten<-Table$Kosten/Table$N
  Table$Zuw.MRSA<-Table$Zuw.MRSA/Table$N
  Table$Zuw.Reg<-Table$Zuw.Reg/Table$N
  Table$Dif<-Table$Zuw.Reg-Table$Zuw.MRSA
  Table$Deck.MRSA<-Table$Zuw.MRSA/Table$Kosten
  Table$Deck.Reg<-Table$Zuw.Reg/Table$Kosten
  Table$Versicherte<-Table$N
  Table$N<-NULL

  writeData(wb, "HMG", Table, rowNames = TRUE)

  ##### Region - wie Oben

  data("OSM_Kreis")           # Lade die Referenztabelle plz zu Kreis
  RS<-data.frame(OSM_Kreis$RS[match(morbiSet@plz,OSM_Kreis$plz)])  # Matching von plz zu Kreis und in Objekt data.frame kopieren
  names(RS)<-"RS"             # Nenne Kreis RS

  Kreis<-sparse.model.matrix(~factor(RS)-1,data=RS)
  IN<-morbiSet@plz %in% OSM_Kreis$plz

  addWorksheet(wb, "Kreis")

  Table<-data.frame(Kosten=(LA[IN]%*%Kreis)@x,Zuw.MRSA=(Z.MRSA[IN]%*%Kreis)@x,Zuw.Reg=(Z.Reg[IN]%*%Kreis)@x,N=colSums(Kreis))
  Table$Kosten<-Table$Kosten/Table$N
  Table$Zuw.MRSA<-Table$Zuw.MRSA/Table$N
  Table$Zuw.Reg<-Table$Zuw.Reg/Table$N
  Table$Dif<-Table$Zuw.Reg-Table$Zuw.MRSA
  Table$Deck.MRSA<-Table$Zuw.MRSA/Table$Kosten
  Table$Deck.Reg<-Table$Zuw.Reg/Table$Kosten
  Table$Versicherte<-Table$N
  Table$N<-NULL

  data("Geo_Ref")           # Lade Geo_ref (Datei mit Namen der Kreise)

  rownames(Table)<-substring(colnames(Kreis),11,25)                      # Name ist Kreisnummer / substring da Name aktuell "factor(RS)16076"
  Table$Region<-Geo_Ref$Raumeinheit[match(rownames(Table),Geo_Ref$RS)]   # Packe hinten an die Tabelle die Namen nach Regionalschlüssel join

  writeData(wb, "Kreis", Table, rowNames = TRUE)



  saveWorkbook(wb, file, overwrite = TRUE)    # Schreibe Excel auf die Festplatte und wenn schon vorhanden dann überschreiben

}
