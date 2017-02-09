\name{SSR}
\alias{SSR}
\title{Shift-Share-Regression}
\usage{
SSR(Data,D,Region.Gewicht="Anteile"\cr,hpruef=FALSE,COV=TRUE,FTest=FALSE\cr,Hochkosten="kein",...)

Data          : Daten vom Typ morbiSet.\cr
D             : (nXk) Matrix mit k Kategorialvariablen.
                Jede Kategorie erfuellt die Nebenbedingung
                ueber alle Beobachtungen.\cr
modell        : Vorimplementiertes Modell fuer D
                (Angabe von D entfaellt).
                Zur Auswahl stehen:
                Kreis: Landkreisdummys.
                Siedlung: Siedlungstyp nach BBSR.
                Raumordnung: Raumordnungsregion nach BBSR.
                Stadt_Land: 1 / 0 fuer Stadt und Land.
                dutch: Regionaldummys des Niedreländischem RSA.\cr
Region.Gewicht: Legt fest ob Nebenbedingung Summarisch
                oder Relativ (Anteile) erfuellt wird.\cr
COV            : Schaetze Varianz-Kovarianz-Matrix.\cr
FTest          : Schaetze einen F-Test.\cr
Hochkosten     : Seperate Hochkostenbewertung
                 (Kein, RP1, RP2, HCP, KP100)

}
\description{
Durchfuehren einer Morbi-RSA Regression mittels Shift-Share-Regression (Gewichtete kleinste Quadrateverfahren).
Schaetzverfahren nach  \link[https://www.jstor.org/stable/2109587]{Greene und Seaks (1991)}

\deqn{Kost = X\beta+\sum Z\gamma+U} mit
\deqn{\sum \gamma = 0}
\deqn{U~(0,1\tage)}

Aktuell ist keine Bewertung von negativen Regressionskoeffizienten implementiert.

Data muss vom Typ \link[MorbiSet]{morbiSet} sein.

}
\examples{
Data<-gen_sample()
# Zufaellige Raumeinteilung (sample)
D<-as.matrix(as.factor(sample(1:100,length(Data@Kost),replace=TRUE))
fit<-SSR(Data,D)

summary(fit)

######## Mit Vorimplementierung

fit.Kreis<-SSR(Data,modell="Kreis")
fit.dutch<-SSR(Data,modell="dutch")

## Excel mit Vergleichszahlen
write.fit(Data,fit.Kreis,fit.dutch,file="Vergleich.xlsx")

}
