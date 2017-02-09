\name{HK}
\alias{HK}
\title{Hochkosten}
\usage{
RP1(Data,HK_Grenze, HK_Ausgleich)

Data          : Daten vom Typ morbiSet.\cr
HK_Grenze     : Untere Grenze fuer Hochkostenfaelle\cr
HK_Ausgleich  : Ausgleichsfaktor
}
\description{
Hochkostenmodelle sind Linkfunktionen, die mit allen Regressionsfunktionen im Paket kombinierbar sind. Innerhalb der Regression ist z.B. "Hochkostenmodell=RP1" anzugeben. Die Funktionen haben ausserhalb der Regression nur geringen Nutzen.\cr

RP1  : Risikopool nach Wasem / Schillo (Standardwerte sind Grenze von 30.000 und Ausgleich von 0,8)
RP2  : Huber-Gewicht (Standardwerte sind Grenze von 1e+07)
HCP  : Capitation-Modell (keine Standardwerte nötig, aber möglich)
KP100 : Trunkierung (Standardwert Grenze von 100.000)

}
\examples{
Data<-gen_sample()
fit<-MRSA(Data,Hochkosten="RP1",HK_Grenze=100000,HK_Ausgleich=0.6)

summary(fit)
}
