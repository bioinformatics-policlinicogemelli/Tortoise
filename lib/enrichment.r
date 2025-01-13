library(enrichR)
options(enrichR.base.address="https://maayanlab.cloud/Enrichr/")
options(enrichR.quiet=TRUE)

all_analisi<- function(params){
  file_geni <- params[[1]]
  cluster <- params[[2]]
  out_path <- params[[3]]

  temp_read<-read.table(
    file_geni,      # TXT data file indicated as string or full path to the file
    header = FALSE, # Whether to display the header (TRUE) or not (FALSE)
    sep = "",       # Separator of the columns of the file
  )
  lista_geni<-list()
  for (l in temp_read){
    lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
  }
  #inseriamo i GENI per l'analisi GO
  response_analisi<- enrichr(genes = lista_geni, databases = c("GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023","KEGG_2021_Human","WikiPathway_2023_Human","Reactome_2022"))

  if (length(response_analisi$GO_Biological_Process_2023) == 1) {
    #print("API not reponde, retrying")
    all_analisi(params)
  } else {
    write.csv(response_analisi$GO_Biological_Process_2023, file = paste(out_path,"/GO/biological_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$GO_Molecular_Function_2023,file = paste(out_path,"/GO/molecular_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$GO_Cellular_Component_2023,file = paste(out_path,"/GO/cellular_",cluster,".csv",sep=""),row.names = FALSE)
    write.csv(response_analisi$KEGG_2021_Human, file = paste(out_path,"/KEGG/kegg_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$WikiPathway_2023_Human, file = paste(out_path,"/WIKI/wiki_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$Reactome_2022, file = paste(out_path,"/REACTOME/reactome_",cluster,".csv",sep=""), row.names = FALSE)
  }
}

#Example analisi_enrich_GO("/Gene/genes_cluster_1.csv", 1, "./output/lung")
analisi_enrich_GO<- function(params){
  file_geni <- params[[1]]
  cluster <- params[[2]]
  out_path <- params[[3]]

  temp_read<-read.table(
    file_geni,      # TXT data file indicated as string or full path to the file
    header = FALSE, # Whether to display the header (TRUE) or not (FALSE)
    sep = "",       # Separator of the columns of the file
  )
  lista_geni<-list()
  for (l in temp_read){
    lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
  }
  #inseriamo i GENI per l'analisi GO
  upEnriched_go<- enrichr(genes = lista_geni, databases = c("GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023"))

  #dati_Biological_Process
  g_process=upEnriched_go[[1]]
  #dati_Molecular
  g_molecular=upEnriched_go[[2]]
  #dati_Cellular
  g_cellular=upEnriched_go[[3]]
  #dataset GO_PROCESS
  write.csv(g_process, file = paste(out_path,"/biological_",cluster,".csv",sep=""), row.names = FALSE)
  #dataset GO_MOLECULAR
  write.csv(g_molecular,file =paste(out_path,"/molecular_",cluster,".csv",sep=""), row.names = FALSE)
  #dataset GO_CELLULAR
  write.csv(g_cellular,file =paste(out_path,"/cellular_",cluster,".csv",sep=""),row.names = FALSE)
}

#Example analisi_enrich_KEGG("/Gene/genes_cluster_1.csv", 1, "./output/lung")
analisi_enrich_KEGG <- function(params) {
  file_geni <- params[[1]]
  cluster <- params[[2]]
  out_path <- params[[3]]

  temp_read<-read.table(
    file_geni,         # TXT data file indicated as string or full path to the file
    header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
    sep = "",          # Separator of the columns of the file
  )
  lista_geni<-list()
  for (l in temp_read){
    lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
  }
  #inseriamo i GENI per l'analisi KEGG
  upEnriched_kegg<- enrichr(genes = lista_geni, databases = c("KEGG_2021_Human"))
  if (length(upEnriched_kegg$KEGG_2021_Human) == 1) {
    print("API not reponde, retrying")
    analisi_enrich_KEGG(params)
  } else {
    #dataset KEGG
    write.csv(upEnriched_kegg, file = paste(out_path,"/kegg_",cluster,".csv",sep=""), row.names = FALSE)
  }
}

#Example analisi_enrich_wiki("/Gene/genes_cluster_1.csv", 1, "./output/lung")
analisi_enrich_wiki<- function(params){
  file_geni <- params[[1]]
  cluster <- params[[2]]
  out_path <- params[[3]]

  temp_read<-read.table(
    file_geni,         # TXT data file indicated as string or full path to the file
    header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
    sep = "",          # Separator of the columns of the file
  )
  lista_geni<-list()
  for (l in temp_read){
    lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
  }
  #inseriamo i GENI per l'analisi WIKI
  upEnriched_wiki<- enrichr(genes = lista_geni, databases =c("WikiPathway_2023_Human") )
  if (length(upEnriched_wiki$WikiPathway_2023_Human) == 1) {
    print("API not reponde, retrying")
    analisi_enrich_wiki(params)
  } else {
    #dataset WIKI
    write.csv(upEnriched_wiki, file = paste(out_path,"/wiki_",cluster,".csv",sep=""), row.names = FALSE)
  }
}

#Example analisi_enrich_reactome("/Gene/genes_cluster_1.csv", 1, "./output/lung")
analisi_enrich_reactome<- function(params){
  file_geni <- params[[1]]
  cluster <- params[[2]]
  out_path <- params[[3]]

  temp_read<-read.table(
    file_geni,         # TXT data file indicated as string or full path to the file
    header = FALSE,    # Whether to display the header (TRUE) or not (FALSE)
    sep = "",          # Separator of the columns of the file
  )  
  lista_geni<-list()
  for (l in temp_read){
    lista_geni<-c(lista_geni,unlist(strsplit(l,",")))
  }
  #inseriamo i GENI per l'analisi REACTOME
  upEnriched_reactome<- enrichr(genes = lista_geni, databases =c("Reactome_2022"))
  if (length(upEnriched_reactome$Reactome_2022) == 1) {
    print("API not reponde, retrying")
    analisi_enrich_reactome(params)
  } else {
    #dataset REACTOME
    write.csv(upEnriched_reactome, file = paste(out_path,"/reactome_",cluster,".csv",sep=""), row.names = FALSE)
  }
}
