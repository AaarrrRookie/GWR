#########################################################################

#                         Methoden für fit

#########################################################################

# Zuweisung der Informationen aus dem Objekt morbiSet
setGeneric(name="setOld",
           def=function(theObject,theOld,name)
           {
             standardGeneric("setOld")
           }
)

setMethod(f="setOld",
          signature="fitMorbiRSA",
          definition=function(theObject,theOld,name)
          {
            theObject@Modell<-name               # Name des Modells setzen
            theObject@Kost <- theOld@Kost        # Setze Kosten
            theObject@tage <- theOld@tage        # Setze Tage
            theObject@plz <- theOld@plz          # Setze Postleitzahl

            return(theObject)                    # Rückgabe an das Objekt
          }
)



#############################

# Zuweisung aus der Regression setzen
setGeneric(name="setZuweisung",
           def=function(theObject,Zuweisung)
           {
             standardGeneric("setZuweisung")
           }
)

setMethod(f="setZuweisung",
          signature="fitMorbiRSA",
          definition=function(theObject,Zuweisung)
          {
            theObject@Zuweisung<-Zuweisung                     # Setze Zuweisung

            R2<-var(Zuweisung)/var(theObject@Kost)             # Berechne R2 usw.
            MAPE<-sum(abs(theObject@Kost-Zuweisung))/length(theObject@Kost)
            CPM<-1-sum(abs(theObject@Kost-Zuweisung))/sum(abs(theObject@Kost-mean(theObject@Kost)))
            Guete<-data.frame(R2,MAPE,CPM)                     # Packe alles in data.frame Guete

            theObject@Guete<-Guete

            return(theObject)                                  # Zurück an das Objekt
          }
)

#############################

# Setze Regressionskoeffizienten im Objekt
setGeneric(name="setCoef",
           def=function(theObject,V,B)
           {
             standardGeneric("setCoef")
           }
)

setMethod(f="setCoef",
          signature="fitMorbiRSA",
          definition=function(theObject,V,B)
          {
            Coef<-data.frame(Coef=names(B),B=B,RF=B/mean(theObject@Kost),SD=sqrt(V),t=B/sqrt(V))    # Sezte Koeffizienten, Risikofaktoren, Standardabweichung und t-Test
            n<-length(theObject@Kost)
            if(n==10){n<-Inf}
            Coef$p<-1-pt(abs(Coef$t),n-length(B))                                   # p-Wert für t-Test
            Coef$star<-""                                                                           # Signifikanz-"Sternchen" für t-Test
            for(i in 1:nrow(Coef)){
              if(!is.na(Coef$p[i])){
                if(Coef$p[i]<0.01){
                  Coef$star[i]<-"***"
                }else if(Coef$p[i]<0.025){
                  Coef$star[i]<-"**"
                }else if(Coef$p[i]<0.05){
                  Coef$star[i]<-"*"
                }
              }
            }
            names(Coef)[1]<-"coef"                                                                  # Setze Namen der Koeffizienten

            theObject@Coef<-Coef                                                                    # Binde an das Objekt

            return(theObject)                                                                       # Rückgabe des Objekt
          }
)


#############################

# Setze F-Test im Objekt
setGeneric(name="setTest",
           def=function(theObject,Z0,Z,B,Dim.X)
           {
             standardGeneric("setTest")
           }
)

setMethod(f="setTest",
          signature="fitMorbiRSA",
          definition=function(theObject,Z0,Z,B,Dim.X)
          {
            F<-((sum((theObject@Kost-Z0)^2)-sum((theObject@Kost-Z)^2))/(length(B)-Dim.X[2]))/(sum((Kost-Z)^2)/freiheitsG)    # F-Test-Statistik siehe Standardlehrbuch
            Ftest<-1-pf(F,df1=(length(B)-Dim.X[2]),df2=length(theObject@Kost)-length(B))                                     # Durchführen des F-Test (p-Wert)


            theObject@FTest<-data.frame(F=F,p=Ftest)                                                                         # Rückgabe des F-Test

            return(theObject)                                                                                                # Rückgabe des Objekt
          }
)



