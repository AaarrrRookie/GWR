\name{MRSA}
\alias{Morbi-RSA}
\title{Morbi-RSA-Regression}
\usage{
MRSA(Data,Standfehler, Hochkosten)

Data          : Daten vom Typ morbiSet.\cr
Standfehler   : Schaetze Varianz-Kovarianz-Matrix.\cr
Hochkosten     : Seperate Hochkostenbewertung
                 (Kein, RP1, RP2, HCP, KP100)
}
\description{
Durchfuehren einer Morbi-RSA Regression (Gewichtete kleinste Quadrateverfahren).

\deqn{Kost = X\beta+U}
\deqn{U~(0,\frac{1}{tage})}

Aktuell ist keine Bewertung von negativen Regressionskoeffizienten implementiert.

Data muss vom Typ \link[MorbiSet]{morbiSet} sein.

Sandfehler  Soll der Standardfehler mit ausgewertet werden (default = TRUE)

}
\examples{
Data<-gen_sample()
fit<-MRSA(Data)

summary(fit)
}
