\name{GWR}
\alias{GWR-RSA}
\title{Geographisch Gewichtete Morbi-RSA-Regression}
\usage{
GWR(Data,k_Krit,...)

Data          : Daten vom Typ morbiSet.
Hochkosten    : Seperate Hochkostenbewertung
                (Kein, RP1, RP2, HCP, KP100)

}
\description{
Durchfuehren einer Morbi-RSA Regression mit geographisch gewichtetem Kerndichteschätzer.
Schaetzverfahren nach  \link[http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471496162.html]{A. Stewart Fotheringham, Chris Brunsdon, Martin Charlton (2002)}

\deqn{Kost = X\beta+\gamma+U}
\deqn{U~(0,\frac{1}{tage})}

Aktuell ist keine Bewertung von negativen Regressionskoeffizienten implementiert.

Data muss vom Typ \link[MorbiSet]{morbiSet} sein.

}
\examples{
Data<-gen_sample()
fit<-GWR(Data,k_Krit=1000)

summary(fit)
}
