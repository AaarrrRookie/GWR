# Definition der Funktion für ein fit
# Funktion entspricht dem Morbi-RSA mit GWR-Regionalfaktor
# angewendet auf das fit nach Kost = Xb + GWR + Fehler

GWR_Deck<-function(fit,k_Krit,std_exakt=FALSE,Kernel="GWR",...){	# GWR-Regression

  # fit - Daten der classe fit
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

  IN<-fit@plz %in% OSM$plz				# Nur Versicherte mit bekanntem Wohnort verwenden, alle Anderen als Wahr in dem Vektor IN markieren

  plz<-as.integer(fit@plz[IN])   # Kopiere in plz von denjenigen die bekannt sind -< als integer

  print("Distanz")

  Dist<-distanz(OSM,fit@plz[IN])			          # Distanz aller Versicherten zueinander mittels distanzfunktion

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

  Zuweisung<-fit@Zuweisung
  Kosten<-fit@Kost

  coef.REG<-GWR_core1_Deck(as.integer(Dist$Match-1),W,Zuweisung,Kosten)
  Erg<-data.frame(plz=Dist$Dplz,Wert=coef.REG)

  return(Erg)
}
