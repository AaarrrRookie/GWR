gen_data_to_sample<-function(morbiSet){
  Data<-data.frame(as.matrix(morbiSet@X))

  Names<-names(Data)
  AGG<-which(substring(Names,1,3)=="AGG"|substring(Names,1,3)=="agg")
  number.AGG<-substring(Names[AGG],4,6)
  number.AGG[nchar(number.AGG)<2]<-paste0("0",number.AGG[nchar(number.AGG)<2])
  Names[AGG]<-paste0("AGG",number.AGG)

  HMG<-which(substring(Names,1,3)=="HMG"|substring(Names,1,3)=="hmg")
  number.HMG<-substring(Names[HMG],4,7)
  for(i in 1:2){
    number.HMG[nchar(number.HMG)<3]<-paste0("0",number.HMG[nchar(number.HMG)<3])
  }
  Names[HMG]<-paste0("HMG",number.HMG)

  EMR<-which(substring(Names,1,3)=="EMR"|substring(Names,1,3)=="emr"|substring(Names,1,3)=="EMG"|substring(Names,1,3)=="emg")
  number.HMG<-substring(Names[EMR],4,7)
  for(i in 1:2){
    number.HMG[nchar(number.HMG)<3]<-paste0("0",number.HMG[nchar(number.HMG)<3])
  }
  Names[EMR]<-paste0("HMG",number.HMG)

  KEG<-which(substring(Names,1,3)=="KEG"|substring(Names,1,3)=="keg")
  number.HMG<-substring(Names[EMR],4,7)
  for(i in 1:2){
    number.HMG[nchar(number.HMG)<3]<-paste0("0",number.HMG[nchar(number.HMG)<3])
  }
  Names[KEG]<-paste0("KEG",number.HMG)

  names(Data)<-Names

  Data$verstorben<-morbiSet@verstorben
  nHMG<-rowSums(Data[,HMG])

  for(i in 0:25){
    Name<-paste0("NHMG",i)
    x<-as.numeric(nHMG==i)
    #attach(Data)
    #assign(Name, x)
    Data[[Name]]<-x
  }

  data("OSM_Kreis")
  M<-match(morbiSet@plz,OSM_Kreis$plz)
  PLZ<-data.frame(OSM_Kreis$RS[M])
  names(PLZ)<-"REG"
  PLZ$REG<-as.character(PLZ$REG)
  PLZ$REG[nchar(PLZ$REG)==4]<-paste0("0",PLZ$REG[nchar(PLZ$REG)==4])
  PLZ$REG[is.na(PLZ$REG)]<-paste0("REG","00000")

  Z<-data.frame(model.matrix(~REG-1,data=PLZ))
  Data<-cbind(Data,Z)

  return(Data)

}
