\name{MoranI}
\alias{MoransI}
\title{MoransI}
\usage{
MoranI.Test(Data,fit,Kernel,k_krit,type,n_sim)

Data          : Daten vom Typ morbiSet.
fit           : Ergebnisobjekt nach Regression (fit)


Kernel        : Wahl der Distanzgewichtung (bi,inf,GWR)\cr
krit          : kritischer Wert\cr
type          : Annahme an die Verteilung des Morans'I. Bei Normalverteilung (normal)
                folgt die Berechnung der Standardabwichung und des standardisierten
                Morans'I der Literatur. Ohne Verteilungsannahme (MC) wird eine
                Monte Carlo Simulation durchgefuehrt (dauert lang) und die Verteilung
                aproximiert. Standardmaeßig ist "normal",
                die Normalverteilungsannahme eingetagen.\cr
n_sim         : Anzahl der Simulationen fuer MC\cr


MoranI(Match,plz,Dplz,D,Y,TSS)

Match         : Zuordnungsvektor zwischen Dplz und Zielvariable
                (z.B. Versicherter 1 hat die 5te Dplz)\cr
plz           : Vektor Postleitzahlen der Versicherten\cr
Dplz          : Vektor der Postleitzahlen in Deutschland\cr
D             : Distanz oder Gewichtungsmatrix der
                Dimension Laenge Dplz X Laenge Dplz\cr
Y             : Vektor der Zielvariable\cr
TSS           : Total sum of squares von Y

}
\description{

Morans'I berechnet die geographische Autokorrelation eines Vektors. Der MI ist definiert als:

\deqn{I = E' W E} wobei die Matrix W durch den Kernel bestimmt wird. Aktuell stehen zur Auswahl:
\deqn{bi: 1-(1-d^2 / krit)^2}
\deqn{inf: 1/d^2}
\deqn{GWR: 1-(1-d^2 / k^2)^2}
wobei GWR über die k-nächsten Nachbarn und keine kritsche Distanz definiert wird.

}
\examples{
Data<-gen_sample()

fit<-MRSA(Data)
MoranI.Test(Data,fit,Kernel="bi",krit=5000)

}
