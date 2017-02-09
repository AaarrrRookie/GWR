coef.fitMorbiRSA<-function(theObject){
  Coef.fit<-theObject@Coef[,2]/365
  names(Coef.fit)<-theObject@Coef[,1]
  return(Coef.fit)

}


predict.fitMorbiRSA<-function(theObject,morbiSet,modell="NA",Hochkosten="kein",...){

  Data<-morbiSet

  if(Hochkosten=="kein"){
    Z_HK<-rep(0,length(Data@Kost))
  }else if(Hochkosten=="KP100"){
    HK<-KP100(Data,...)
    Z_HK<-HK$HK_Z
  }else if(Hochkosten=="RP1"){
    HK<-RP1(Data,...)
    Z_HK<-HK$HK_Z
  }else if(Hochkosten=="RP2"){
    HK<-RP2(Data,...)
    Z_HK<-HK$HK_Z
  }else if(Hochkosten=="HCP"){
    HK<-HCP(Data,...)
    Z_HK<-HK$HK_Z
  }else{
    stop("Keine Risikopool definiert")
  }





  if(theObject@Modell=="MRSA"){
    Coef.MRSA<-coef(theObject)

    neue_Daten<-morbiSet@X@Dimnames[[2]]
    Matching_Daten<-match(neue_Daten,names(Coef.MRSA))

    if(length(names(Coef.MRSA))>length(neue_Daten)){
      Daten_fehlen<-names(Coef.MRSA)[which(!names(Coef.MRSA) %in% neue_Daten)]
      warning(paste0("Daten Fehlen: ",Daten_fehlen," ;"))

    }

    return((morbiSet@X%*%Coef.MRSA[Matching_Daten])@x*morbiSet@tage+Z_HK)

  }else if(theObject@Modell=="GWR"){
    Coef.GWR<-coef(theObject)

    neue_Daten<-morbiSet@X@Dimnames[[2]]
    Matching_Daten<-match(neue_Daten,names(Coef.GWR))

    Yhat.MRSA<-(morbiSet@X%*%Coef.GWR[Matching_Daten])@x*morbiSet@tage

    plz<-morbiSet@plz
    neue_plz<-unique(plz)

    Matching_Daten<-match(plz,names(Coef.GWR))

    Yhat.GWR<-Coef.GWR[Matching_Daten]*morbiSet@tage



    if(length(names(Coef.GWR))>length(neue_Daten)+length(neue_plz)){
      Daten_fehlen<-c(names(Coef.GWR)[which(!names(Coef.GWR) %in% neue_Daten)],neue_plz[which(!neue_plz %in% names(Coef.GWR))])
      warning(paste0("Daten Fehlen: ",Daten_fehlen," ;"))

    }

    return(Yhat.MRSA+Yhat.GWR+Z_HK)

  }else if(theObject@Modell=="SSR"){

    Data<-morbiSet

    Coef.SSR<-coef(theObject)

    neue_Daten<-morbiSet@X@Dimnames[[2]]
    Matching_Daten<-match(neue_Daten,names(Coef.SSR))

    Yhat.MRSA<-(morbiSet@X%*%Coef.SSR[Matching_Daten])@x*morbiSet@tage

    if(modell=="NA"){

    }else if(modell=="Kreis"){
      data("OSM_Kreis")
      RS<-data.frame(OSM_Kreis$RS[match(Data@plz,OSM_Kreis$plz)])
      names(RS)<-"RS"
      RS[is.na(RS)]<-00000

      D<-as.matrix(as.factor(RS$RS))

    }else if(modell=="Bundesland"){
      data("OSM_Kreis")
      RS<-data.frame(OSM_Kreis$RS[match(Data@plz,OSM_Kreis$plz)])
      names(RS)<-"RS"
      RS$RS<-as.character(RS$RS)
      RS$RS[nchar(RS$RS)==4]<-paste0("0",RS$RS[nchar(RS$RS)==4])
      RS$RS<-as.numeric(substring(RS$RS,1,2))
      RS[is.na(RS)]<-0

      D<-as.matrix(as.factor(RS$RS))

    }else if(modell=="dutch"){
      data("OSM_Kreis")
      data("Geo_Ref")
      RS<-data.frame(OSM_Kreis$RS[match(Data@plz,OSM_Kreis$plz)])
      names(RS)<-"RS"

      RS$RS<-Geo_Ref$dutch[match(RS$RS,Geo_Ref$RS)]
      RS[is.na(RS)]<-00000

      D<-as.matrix(as.factor(RS$RS))
    }else if(modell=="Siedlung"){
      data("OSM_Kreis")
      data("Geo_Ref")
      RS<-data.frame(OSM_Kreis$RS[match(Data@plz,OSM_Kreis$plz)])
      names(RS)<-"RS"

      RS$RS<-Geo_Ref$Siedlung[match(RS$RS,Geo_Ref$RS)]
      RS[is.na(RS)]<-00000

      D<-as.matrix(as.factor(RS$RS))
    }else if(modell=="Raumordnung"){
      data("OSM_Kreis")
      data("Geo_Ref")
      RS<-data.frame(OSM_Kreis$RS[match(Data@plz,OSM_Kreis$plz)])
      names(RS)<-"RS"

      RS$RS<-Geo_Ref$Raumordnung[match(RS$RS,Geo_Ref$RS)]
      RS[is.na(RS)]<-00000

      D<-as.matrix(as.factor(RS$RS))
    }else if(modell=="Stadt_Land"){
      data("OSM_Kreis")
      data("Geo_Ref")
      RS<-data.frame(OSM_Kreis$RS[match(Data@plz,OSM_Kreis$plz)])
      names(RS)<-"RS"

      RS$RS<-Geo_Ref$Stadt_Land[match(RS$RS,Geo_Ref$RS)]
      RS[is.na(RS)]<-00000

      D<-as.matrix(as.factor(RS$RS))
    }else{
      stop("keine Zusatzinformation")
    }

    for(i in 1:ncol(D)){

      Z<-sparse.model.matrix(~as.factor(D[,i])-1)
    }

    neue_Daten<-Z@Dimnames[[2]]
    Matching_Daten<-match(neue_Daten,names(Coef.SSR))


    Yhat.SSR<-(Z%*%Coef.SSR[Matching_Daten])@x*morbiSet@tage

    return(Yhat.MRSA+Yhat.SSR+Z_HK)

  }


}
