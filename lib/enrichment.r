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
  response_analisi<- enrichr(genes = lista_geni, databases = c("GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023","KEGG_2021_Human","WikiPathway_2023_Human","Reactome_2022"))
  if (length(response_analisi$GO_Biological_Process_2023) == 1) {
    # API not reponde, retrying
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
