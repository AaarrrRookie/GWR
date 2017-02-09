metaReg<-function(list_of_files,sheet){

  object_lenge<-c()
  Namen_liste<-c()

  for(file in list_of_files){
     Name<-substring(file,1,nchar(file)-5)
     Namen_liste<-c(Namen_liste,Name)

     assign(Name,read.xlsx(file,sheet=sheet))
     object_lenge<-c(object_lenge,nrow(get(Name)))
  }

  start<-Namen_liste[which(object_lenge==max(object_lenge))[1]]

  Table<-get(start)
  Table<-data.frame(Table$coef,Table$B,1/Table$SD)
  names(Table)<-c("Coef",paste0("B_",start),paste0("Prez_",start))

  Meta<-Table[,2]*Table[,3]
  Prez<-Table[,3]

  for(i in 1:length(Namen_liste)){
    if(!Namen_liste[i]==start){
      Table2<-get(Namen_liste[i])
      M<-match(Table[,1],Table2[,1])
      Table[[paste0("B_",Namen_liste[i])]]<-Table2$B[M]
      Table[[paste0("Prez_",Namen_liste[i])]]<-1/Table2$SD[M]

      Meta<-Meta+Table2$B/Table2$SD[M]
      Prez<-Prez+1/Table2$SD[M]
    }
  }

  Table$B_Meta<-Meta/Prez
  Table$Prez_Meta<-Prez

  wb <- createWorkbook()
  options("openxlsx.borderStyle" = "thin")                # Tabellen-Definition
  options("openxlsx.borderColour" = "#4F81BD")
  addWorksheet(wb, "Meta_Reg")
  writeData(wb, "Meta_Reg", Table, rowNames = TRUE)

  saveWorkbook(wb, "Meta.xlsx", overwrite = TRUE)

}

meta_to_fit<-function(file,Modell){
  fit<-fitMorbiRSA()
  meta<-read.xlsx(file,sheet=1)
  B<-meta$B_Meta
  names(B)<-meta$Coef
  SD<-1/meta$Prez_Meta

  fit<-setCoef(fit,SD^2,B)

  fit@Modell<-Modell
  return(fit)
}
