#install.packages("devtools")
#if(!require(devtools)) install.packages("devtools",repos = "http://cran.rstudio.com/")
#install.packages("install_github")
#remotes::install_github("wjawaid/enrichR")
#install.packages("jsonlite")

library(jsonlite)
library(devtools)
library(ggplot2)
#importiamo la libreria enrichR
library(enrichR)

listEnrichrSites()
#importiamo i termini di enrichR per l'analisi
dbs <- listEnrichrDbs()
#head(dbs)


#ANALISI GO
#creiamo il dbs per i termini GO
#1-->Pocess
#2-->Molecular
#2-->Cellular
dbs_go <- c("GO_Biological_Process_2023","GO_Molecular_Function_2023", "GO_Cellular_Component_2023")
dbs_kegg <- c("KEGG_2021_Human")
dbs_wiki <- c("WikiPathway_2023_Human")
dbs_reactome <- c("Reactome_2022")

json_data <- fromJSON("./config.json")

PATH=json_data$Paths$output_folder
go=json_data$Enrichment$go
kegg=json_data$Enrichment$kegg
wiki=json_data$Enrichment$wiki
reactome=json_data$Enrichment$reactome

path_gene <- paste0(PATH,"/Gene")
file_list <- list.files(path_gene)
# Conta il numero di file
num_files <- length(file_list)
N_CLUSTER=(num_files-1)
#N_CLUSTER=4
if(go==T){
analisi_enrich_GO<- function(file_geni,cluster){
  print(paste0(PATH,file_geni))
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
  
  if (!dir.exists(paste0(PATH,"/Arricchimento_all_genes/GO"))) {
    # Se la cartella non esiste, la creiamo
    dir.create(paste0(PATH,"/Arricchimento_all_genes/GO"),recursive=TRUE)
    cat("Cartella creata:", (paste0(PATH,"/Arricchimento_all_genes/GO")), "\n")
  } #else {
    #cat("La cartella esiste già:", (paste0(PATH,"/Arricchimento_all_genes/GO")), "\n")
  #}
  
  #dataset GO_PROCESS

  write.csv(g_process, file = paste(PATH,"/Arricchimento_all_genes/GO/biological_",cluster,".csv",sep=""), row.names = FALSE)
  #dataset GO_MOLECULAR
  write.csv(g_molecular,file =paste(PATH,"/Arricchimento_all_genes/GO/molecular_",cluster,".csv",sep=""), row.names = FALSE)
  #dataset GO_CELLULAR
  write.csv(g_cellular,file =paste(PATH,"/Arricchimento_all_genes/GO/cellular_",cluster,".csv",sep=""),row.names = FALSE)
  
}


#analisi_enrich_GO

 for (i in 0:N_CLUSTER){
  print(i)
  analisi_enrich_GO(paste0("/Gene/genes_cluster_",i,".csv"),i)}

}
#*******************************

#KEGG
if(kegg==T){
analisi_enrich_KEGG<- function(file_geni,cluster){
  temp_read<-read.table(paste0(PATH,file_geni),            # TXT data file indicated as string or full path to the file
                      header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
                      sep = "",          # Separator of the columns of the file
)  

lista_geni<-list()
for (l in temp_read){
  lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
}

#inseriamo i GENI per l'analisi KEGG
upEnriched_kegg<- enrichr(genes = lista_geni, databases = dbs_kegg)
  
if (!dir.exists(paste0(PATH,"/Arricchimento_all_genes/KEGG"))) {
  # Se la cartella non esiste, la creiamo
  dir.create(paste0(PATH,"/Arricchimento_all_genes/KEGG"),recursive=TRUE)
  cat("Cartella creata:", (paste0(PATH,"/Arricchimento_all_genes/KEGG")), "\n")
} #else {
  #cat("La cartella esiste già:", (paste0(PATH,"/Arricchimento_all_genes/KEGG")), "\n")
#}

  #dataset KEGG
write.csv(upEnriched_kegg, file = paste(PATH,"/Arricchimento_all_genes/KEGG/kegg_",cluster,".csv",sep=""), row.names = FALSE)
 
}

#analisi_enrich_KEGG
 for (i in 0:N_CLUSTER){
  print(i)
  analisi_enrich_KEGG(paste0("/Gene/genes_cluster_",i,".csv"),i)
 }

}

#************************************

 #WikiPathway
if(wiki==T){
analisi_enrich_wiki<- function(file_geni,cluster){
  temp_read<-read.table(paste0(PATH,file_geni),            # TXT data file indicated as string or full path to the file
                      header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
                      sep = "",          # Separator of the columns of the file
)  

lista_geni<-list()
for (l in temp_read){
  lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
}

#inseriamo i GENI per l'analisi WIKI
upEnriched_wiki<- enrichr(genes = lista_geni, databases =dbs_wiki )


if (!dir.exists(paste0(PATH,"/Arricchimento_all_genes/WIKI"))) {
  # Se la cartella non esiste, la creiamo
  dir.create(paste0(PATH,"/Arricchimento_all_genes/WIKI"),recursive=TRUE)
  cat("Cartella creata:", (paste0(PATH,"/Arricchimento_all_genes/WIKI")), "\n")
} #else {
  #cat("La cartella esiste già:", (paste0(PATH,"/Arricchimento_all_genes/WIKI")), "\n")
#}
  
  #dataset WIKI
write.csv(upEnriched_wiki, file = paste(PATH,"/Arricchimento_all_genes/WIKI/wiki_",cluster,".csv",sep=""), row.names = FALSE)
 
}

#analisi_enrich_WIKI
 for (i in 0:N_CLUSTER){
  print(i)
  analisi_enrich_wiki(paste0("/Gene/genes_cluster_",i,".csv"),i)
}
}
#*****************************+
#reactome
if(reactome==T){
 analisi_enrich_reactome<- function(file_geni,cluster){
  temp_read<-read.table(paste0(PATH,file_geni),            # TXT data file indicated as string or full path to the file
                      header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
                      sep = "",          # Separator of the columns of the file
)  

lista_geni<-list()
for (l in temp_read){
  lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
}

#inseriamo i GENI per l'analisi REACTOME
upEnriched_reactome<- enrichr(genes = lista_geni, databases =dbs_reactome)

if (!dir.exists(paste0(PATH,"/Arricchimento_all_genes/REACTOME"))) {
  # Se la cartella non esiste, la creiamo
  dir.create(paste0(PATH,"/Arricchimento_all_genes/REACTOME"),recursive=TRUE)
  cat("Cartella creata:", (paste0(PATH,"/Arricchimento_all_genes/REACTOME")), "\n")
} #else {
 # cat("La cartella esiste già:", (paste0(PATH,"/Arricchimento_all_genes/REACTOME")), "\n")
#}
  
  #dataset REACTOME
  write.csv(upEnriched_reactome, file = paste(PATH,"/Arricchimento_all_genes/REACTOME/reactome_",cluster,".csv",sep=""), row.names = FALSE)
 
}

#analisi_enrich_REACTOME
 for (i in 0:N_CLUSTER){
  print(i)
   analisi_enrich_reactome(paste0("/Gene/genes_cluster_",i,".csv"),i)
 }
}
 