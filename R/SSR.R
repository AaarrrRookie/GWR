SSR<-function(Data                         # Daten
              ,D=NA                        # Zusatzdaten z.B. Kreis
              ,modell="NA"                 # Zusatzmodell aktuell nur Kreis, Siedlung, Raumordnung, Stadt_Land und dutch möglich
              ,Region.Gewicht="Anteile"    # soll pro Kopf oder pro VT (Anteile) gewichtet werden
              ,hpruef=FALSE                # Hypothesentests (ja / nein)
              ,COV=TRUE                    # Varianz-Covarianz (ja / nein)
              ,FTest=FALSE                 # F-Test (ja / nein)
              ,Hochkosten="kein"           # Hochkostenmodell
              ,...                         # Weitere Einstellungen wenn Hochkostenmodell gewählt
){ # shift and share regression

  # Kost 	= Kosten
  # X		= Designmatrix
  # D		= Matrix mit zusätzlichen (regionalen) Risikofaktoren als Faktor z.B. (Kreis,Bundesland) für Kreis + Bundeslandeffekt in einem Model
  # tage	= Versichertentage
  # Region.Gewicht = "Anteile" für Versichertenage als Regressionsbedingung sonst Köpfe
  # 		Hypotesentest, Covarianz, Ftest

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


  fit<-fitMorbiRSA()
  fit<-setOld(fit,Data,"SSR")

  Dim.X<-dim(Data@X)


	solver<-function(X,Y,w,M){
		Xw<-X
		Xw@x<-w[M]
		return(solve(crossprod(Xw,X),crossprod(Xw,Y)))
	}

	tage<-Data@tage
	w<-tage/max(tage)

	M<-Data@X@i+1

	if(Hochkosten=="kein"){
	  Y<-Data@Kost/w
	  Z_HK<-rep(0,length(Y))
	}else if(Hochkosten=="KP100"){
	  HK<-KP100(Data,...)
	  Y<-HK$LA_NK/w
	  Z_HK<-HK$HK_Z
	}else if(Hochkosten=="RP1"){
	  HK<-RP1(Data,...)
	  Y<-HK$LA_NK/w
	  Z_HK<-HK$HK_Z
	}else if(Hochkosten=="RP2"){
	  HK<-RP2(Data,...)
	  Y<-HK$LA_NK/w
	  Z_HK<-HK$HK_Z
	}else if(Hochkosten=="HCP"){
	  HK<-HCP(Data,...)
	  Y<-HK$LA_NK/w
	  Z_HK<-HK$HK_Z
	}else{
	  stop("Keine Risikopool definiert")
	}




	if(length(w==0)>0){
	  Y[w==0]<-0
	}

	if(FTest==TRUE){

		print("0-Modell")

		B0<-solver(Data@X,Y,w,Data@X@i+1)
		Z0<-Data@X%*%B0*w+Z_HK
	}


#################################################################################################
#                                                                                               #
#                                                                                               #
#                           Nebenbedingung				                                        #
#                                                                                               #
#################################################################################################


	for(i in 1:ncol(D)){

		Z<-sparse.model.matrix(~as.factor(D[,i])-1)

		if(i==1){
			R<-matrix(rep(0,dim(Data@X)[2]),1,dim(Data@X)[2])
			if(Region.Gewicht=="Anteile"){
				r<-matrix((tage%*%Z)/sum(as.numeric(tage)),1,ncol(Z))
			}else{
				r<-matrix(rep(1,ncol(Z)),1,ncol(Z))
			}
			R<-cbind(R,r)
		}else{
			R<-rbind(R,0)
			if(Region.Gewicht=="Anteile"){
				r<-matrix((tage%*%Z)/sum(as.numeric(tage)),1,ncol(Z))
			}else{
				r<-matrix(rep(1,ncol(Z)),1,ncol(Z))
			}
			r<-rbind(0,r)
			R<-cbind(R,r)
		}
		R<<-R


		X<-cBind(Data@X,Z)				# Neue Designmatrix aus X und Z
}

	r<-matrix(0,nrow(R),1)				# Bedingungen sollen 0 sein

	print("Regression")

	### Koeffizienten


	W<-cBind(crossprod(X*w,X),t(R))


	for(i in 1:nrow(R)){
		R<-cbind(R,0)
		}



	W<-rBind(W,R)						# Multiplikatormatrix

	v<-rBind(crossprod(X*w,Y),r)		# Lösungsvektor


#################################################################################################
#                                                                                               #
#                                                                                               #
#                           Regression  				                                        #
#                                                                                               #
#################################################################################################


	vW<-qr.solve(W)						# Inverse mittels QR-Zerlegung
	B<-vW%*%v							# Schätzer + lambda

	B<-B[-length(B)]					# Lagrange brauchen wir nicht

	Coef<-data.frame(B)

  if(COV==FALSE){
    fit<-setCoef(fit,rep(0,length(Data@Kost)),B)
  }

#################################################################################################
#                                                                                               #
#                                                                                               #
#                           Kovarianzmatrix				                                        #
#                                                                                               #
#################################################################################################

if(COV==TRUE){

	print("Varianz-Kovarianz")

	freiheitsG<-nrow(X)-(ncol(X)-nrow(R))

	u<-Y-X%*%B
	sigma<-(t(u)%*%u)/freiheitsG	# Freiheitsgrad Beobachtung n - Variablen + Nebenbedingung
	sigma<-sigma[1,1]


	h1<-matrix(0,nrow(R),ncol(X))

	V<-cBind(sigma*(crossprod(X*w,X)),t(h1))		# mittlere Teil der Varianz-Covarianz-Matrix

	for(i in 1:nrow(R)){
		h1<-cbind(h1,0)
	}

	V<-rBind(V,h1)
	V<-diag(vW%*%V%*%vW)[-(length(B)+1)] 								# vollständige Varianz-Covarianz-Matrix

	names(B)<-colnames(X)
	fit<-setCoef(fit,V,B)

}



	Z<-(X%*%B)@x*w + Z_HK					# Zuweisungsvektor

	fit<-setZuweisung(fit,Z)

	if(abs(sum(Z-Data@Kost))>100){
	  warning(paste0("Fehler in der Regresion: ",sum(Z-Data@Kost)))
	}




	#### F-Test

	if(FTest==TRUE){
		print("F-Test")

	fit<-setTest(fit,Z0,Z,B,dim.X)


	}


return(fit)

}


