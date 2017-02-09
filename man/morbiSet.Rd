\name{morbiSet}
\alias{MorbiSet}
\title{Morbiditaetsorientierte Satzarten}
\usage{
morbiSet(X, Kost, tage, plz)

gen_morbiSet(formula, tage, plz, Data)
}
\description{
Aufbau eines Datensatzes vom Typ "morbiSet". MorbiSets enthalten MRSA Versicherteninformatinen aller Risikogruppen (rechte Handseite) in der Variablen X, Leistungsausgaben als "Kost", die Versichertenepisode und die Postleitzahl.

X       - dgCMatrix (Dummy) enhaelt alle Risikogruppen \cr
Kost    - Leistungsausgaben der HLB 1 - 6 \cr
tage    - Versichertentage des Ausgleichsjahres \cr
plz     - Wohnort des Versicherten \cr

}
\examples{
Data <- morbiSet()


#not run

Data <- data.frame(Kost = rnorm(100,0,1), X = rnorm(100,0,1), plz = 1:100, tage = 365)
formula <- Kost ~ X
Test.morbiSet <- gen_morbiSet(formula,"tage","plz",Data)


}
