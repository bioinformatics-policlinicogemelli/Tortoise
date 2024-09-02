

#install.packages("devtools")
#if(!require(devtools)) install.packages("devtools",repos = "http://cran.rstudio.com/")
library(devtools)
install_github("wjawaid/enrichR")
library(ggplot2)


#importiamo la libreria enrichR
library(enrichR)
listEnrichrSites()
#importiamo i termini di enrichR per l'analisi
dbs <- listEnrichrDbs()
head(dbs)


#ANALISI GO
#creiamo il dbs per i termini GO
#1-->Pocess
#2-->Molecular
#2-->Cellular
dbs_go <- c("GO_Biological_Process_2023","GO_Molecular_Function_2023", "GO_Cellular_Component_2023")
dbs_kegg <- c("KEGG_2021_Human")
dbs_wiki <- c("WikiPathway_2023_Human")
dbs_phengeni <- c("PhenGenI_Association_2021")

PATH = "./../../NEW_DATA_NETWORK/Output/LUNG/Lung_Gene_Particular/EGFR_all/"


name="egfr_del_19"


analisi_enrich_GO<- function(file_geni){
  temp_read<-read.table(paste0(PATH,file_geni),            # TXT data file indicated as string or full path to the file
                      header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
                      sep = "",          # Separator of the columns of the file
)  

lista_geni<-list()
for (l in temp_read){
  lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
}



  
  #inseriamo i GENI per l'analisi GO
  upEnriched_go<- enrichr(genes = lista_geni, databases = dbs_go)
  
  #dati_Biological_Process
  g_process=upEnriched_go[[1]]
  #dati_Molecular
  g_molecular=upEnriched_go[[2]]
  #dati_Cellular
  g_cellular=upEnriched_go[[3]]
  
  #dataset GO_PROCESS

  write.csv(g_process, file = paste(PATH,"Arricchimento_",name,"/","biological.csv",sep=""), row.names = FALSE)
  #dataset GO_MOLECULAR
  write.csv(g_molecular,file =paste(PATH,"Arricchimento_",name,"/","molecular.csv",sep=""), row.names = FALSE)
  #dataset GO_CELLULAR
  write.csv(g_cellular,file =paste(PATH,"Arricchimento_",name,"/","cellular.csv",sep=""),row.names = FALSE)
  
}


#analisi_enrich_GO

analisi_enrich_GO(paste0("genes_",name,".txt"))

#*******************************

