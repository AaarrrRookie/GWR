#####################################################################

#                           Definition und Funktionen für das Objekt fit (Ergebnis einer Regression)

######################################################################
# Erstelle Ergebnisse Fitting
#
# Aufbau des fit für geographische Morbi-RSA
fitMorbiRSA <- setClass(
  # Name der Klasse
  "fitMorbiRSA",

  # Slots der Klasse
  slots = c(
    Modell="character",                # Name des Modells
    Kost = "numeric",                  # Kosten der Versicherten im morbiSet
    Zuweisung = "numeric",             # Zuweisungen der Versicherten
    tage="numeric",                    # Versichertentage
    plz="integer",                     # Postleitzahl der Versicherten
    Coef="data.frame",                 # Koeffizienten der Regression
    Guete="data.frame",                # Guete der Regression
    FTest="data.frame"                 # F-Test falls vorhanden
  ),

  # Standardwerte für Slots
  prototype=list(
    Kost = rep(1,10),
    tage = rep(365,10),
    plz = 1:10
  ),

  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  {
    if(abs(sum(object@Kost) - sum(object@Zuweisung))<1000) {          # Teste of Summe der Kosten gleich Summe der Zuweisung
      return("Kosten ungleich Zuweisung")
    }
    return(TRUE)
  }
)


