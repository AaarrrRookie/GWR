library(Matrix)

######################################################################
# Create the MorbiSet
#
# Aufbau des MorbiSet für geographische Morbi-RSA


morbiSet <- setClass(
  # Name der Klasse = morbiSet
  "morbiSet",

  # Definiere slots der Klasse
  slots = c(
    X = "dgCMatrix",        # dünn besetzte Matrix
    Kost = "numeric",
    tage="numeric",
    verstorben="numeric",
    plz="integer"
  ),

  # Standardwerte der slots
  prototype=list(
    X = as(Matrix(matrix(c(0,1),10,10)),"dgCMatrix"),      # 10x10 Matrix aus 0 und 1
    Kost = rep(1,10),                                      # 1x10 Vektor bestehend aus der Zahl 1
    tage = rep(365,10),# 1x10 Vektor bestehend aus der Zahl 365
    verstorben = rep(1,10),
    plz = 1:10                                             # 1x10 Vektor bestehend aus den Zahlen 1 bis 10
  ),

  # Funtion zur kontrolle der Validität der eingebundenen Daten
  # Nicht ausgeführt wenn eine Initialisierungsfunktion ausgeführt wird!
  validity=function(object)
  {
    if(!length(object@Kost) == nrow(object@X) ||
       !length(object@Kost) == length(object@plz) ||
       !length(object@Kost) == length(object@tage)
       ) {
      return("Objektlaenge ungleich")
    }# Wenn die Objekte unterschiedliche Dimensionen haben wird ein Fehler ausgegeben
    if(det(crossprod(object@X)) == 0){
      return("X nicht invertierbar!")
    }

    return(TRUE)
  }
)


# Definiere eine neue Prozedur für das Objekt Morbiset
setGeneric(name="subsetMorbiSet",    # name der Prozedur ist subsetMorbiset
           def=function(theObject,subset)  # Prozedur ist eine Funktion die das Objekt morbiSet und eine Definition für subset benötigt
           {
             standardGeneric("subsetMorbiSet") # Prozedur wird standardmaeßig für morbiSet angewendet
           }
)

# Definiere eine Methode für die Prozedur subsetMorbiSet
setMethod(f="subsetMorbiSet",   # für was
          signature="morbiSet", # über welches Objekt
          definition=function(theObject,subset) # was wird gemacht: Funktion über Objekt und Subset
          {
            theObject@X <- theObject@X[subset,]      # Ordne dem Objekt nur Elemente aus dem Originalobjekt zu wenn Subset = Falsch für slot X
            theObject@Kost <- theObject@Kost[subset] # für slot Kost
            theObject@tage <- theObject@tage[subset] # für slot tage
            theObject@plz <- theObject@plz[subset] # für slot tage

            return(theObject)
          }
)

# Definiere eine neue Prozedur für das Objekt Morbiset
setGeneric(name="sampleMorbiSet",    # name der Prozedur ist subsetMorbiset
           def=function(theObject,sample)  # Prozedur ist eine Funktion die das Objekt morbiSet und eine Definition für subset benötigt
           {
             standardGeneric("sampleMorbiSet") # Prozedur wird standardmaeßig für morbiSet angewendet
           }
)

# Definiere eine Methode für die Prozedur subsetMorbiSet
setMethod(f="sampleMorbiSet",   # für was
          signature="morbiSet", # über welches Objekt
          definition=function(theObject,sample) # was wird gemacht: Funktion über Objekt und Subset
          {
            theObject@X <- theObject@X[sample,]      # Ordne dem Objekt nur Elemente aus dem Originalobjekt zu wenn Subset = Falsch für slot X
            theObject@Kost <- theObject@Kost[sample] # für slot Kost
            theObject@tage <- theObject@tage[sample] # für slot tage
            theObject@plz <- theObject@plz[sample] # für slot tage

            return(theObject)
          }
)


##################################################################################################################################################



# Initialisierungsfunktion für Objekt morbiSet
gen_morbiSet <- function(formula,tage,plz,verstorben,data){ # Funktion erwartet eine Formel, tage als string, plz als string und eine Herkunft (data) kann SQL, csv oder Rdata sein


  X <- sparse.model.matrix(formula,data)    # X ist dünn besetze Matrix definiert aus der Formel z.B. ~ 1 ergibt eine Matrix mit einer Spalte die nur 1 enthält oder
                                            # ~ hmg1 + hmg2 ergibt eine MAtrix X die aus zwei Spalten mit den Informationen hmg1=1 und hmg2=1 besteht.
  Kost <- data[,as.character(formula[[2]])] # Zuordnung was Kosten sind z.B. hlb1~hmg1 liefert für die linke Handseite hlb1
  tage <- data[,tage]                       # Zuordnung der Tage
  plz  <- data[,plz]                        # Zuordnung der Postleitzahl
  verstorben  <- data[,verstorben]                        # Zuordnung der verstorben

  Names<-X@Dimnames[[2]]
  Names<-gsub('factor(','',Names,fixed=TRUE)
  Names<-gsub(')','',Names,fixed=TRUE)
  Names<-gsub('_','',Names,fixed=TRUE)
  X@Dimnames[[2]]<-Names


  theObject<-morbiSet(X=X,Kost=Kost,tage=tage,plz=plz,verstorben=verstorben) # das Objekt wird als morbiSet mit den Komponenten gebaut

  return(theObject)                        # gebe das Objekt zurück (klasse morbiSet)


}



