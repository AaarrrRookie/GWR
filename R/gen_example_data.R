######################################################################################################
#
#                                     Erstellen von Beispieldaten
#
#####################################################################################################

gen_sample<-function(
  n=100000,   # Anzahl der Versicherten
  k=100      # Anzahl HMG
){

data(OSM)                                     # Läd OSM (PLZ zu Punkten) Referenzdatei

plz<-sample(OSM$plz,n,replace=TRUE)           # Zufallsauswalh an PLZ
n_reg<-length(unique(plz))                    # Wieviele Regionen wurden Ausgewählt (einzigartige PLZ)


# Generiert Designmatrix X aus AGG, HMG und Region als Dummies
X<-sparse.model.matrix(~AGG+HMG+Reg-1,data=data.frame(AGG=as.factor(sample(1:40,n,replace=TRUE)),
                                              HMG=as.factor(sample(0:k,n,replace=TRUE)),
                                              Reg=as.factor(plz)))

# Generiert Zufallsvektor für Kosten pro Tag
Y<-(X%*%c((1:40)*1000,(1:k)*10000,(1:(n_reg-1))*100)+rnorm(n,2500,5000))@x/365

# Generiert Designmatrix X ohne Postleitzahlen
X<-sparse.model.matrix(~AGG+HMG-1,data=data.frame(AGG=as.factor(sample(1:40,n,replace=TRUE)),
                                                      HMG=as.factor(sample(0:k,n,replace=TRUE))))

X<-as(X,"dgCMatrix")

# Zufallsauswahl der Versichertentage
tage<-sample(1:365,n,replace=TRUE)
Y<-Y*tage                                                # Kosten je Versicherten

return(morbiSet(X=X,Kost=Y,tage=tage,plz=plz))           # Rückgabe der Werte
}
