\name{write.fit}
\alias{write.fit}
\title{Schreiben der Ergebnisse}
\usage{
write.fit(morbiSet,fit.MRSA,fit.Reg,file="")


}
\description{
Schreiben der Ergebnisse in eine Excel-Datei
}
\examples{


######## Ergebnisvergleich

fit.Kreis<-SSR(Data,modell="Kreis")
fit.dutch<-SSR(Data,modell="dutch")

## Excel mit Vergleichszahlen
write.fit(Data,fit.Kreis,fit.dutch,file="Vergleich.xlsx")

######## Nur MRSA

fit.MRSA<-MRSA(Data)
write.fit(Data,fit.MRSA=fit.MRSA,file="nur_MRSA.xlsx")

######## Nur Regionalmodell

fit.Kreis<-SSR(Data,modell="Kreis")
write.fit(Data,fit.Reg=fit.Kreis,file="nur_Region.xlsx")

}
