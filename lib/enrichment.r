library(enrichR)
options(enrichR.base.address="http://maayanlab.cloud/Enrichr/")
options(enrichR.quiet=TRUE)

all_analisi<- function(params){
  file_geni <- params[[1]]
  cluster <- params[[2]]
  out_path <- params[[3]]
  event_type <- params[[4]]  # MUT, CNV_GAIN, or CNV_LOSS
  
  temp_read<-read.table(
    file_geni,      # TXT data file indicated as string or full path to the file
    header = FALSE, # Whether to display the header (TRUE) or not (FALSE)
    sep = "\t",     # Separator of the columns of the file
  )
  lista_geni<-as.character(unlist(temp_read))
  
  # Skip if no genes
  if (length(lista_geni) == 0) {
    return()
  }
  
  response_analisi<- enrichr(genes = lista_geni, databases = c("GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023","KEGG_2021_Human","WikiPathway_2023_Human","Reactome_2022"))
  if (length(response_analisi$GO_Biological_Process_2023) < 9) {
    # API not responded, retrying
    all_analisi(params)
  } else {
    # Create event_type subdirectories
    dir.create(paste(out_path,"/",event_type,"/GO",sep=""), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste(out_path,"/",event_type,"/KEGG",sep=""), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste(out_path,"/",event_type,"/WIKI",sep=""), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste(out_path,"/",event_type,"/REACTOME",sep=""), showWarnings = FALSE, recursive = TRUE)
    
    write.csv(response_analisi$GO_Biological_Process_2023, file = paste(out_path,"/",event_type,"/GO/biological_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$GO_Molecular_Function_2023,file = paste(out_path,"/",event_type,"/GO/molecular_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$GO_Cellular_Component_2023,file = paste(out_path,"/",event_type,"/GO/cellular_",cluster,".csv",sep=""),row.names = FALSE)
    write.csv(response_analisi$KEGG_2021_Human, file = paste(out_path,"/",event_type,"/KEGG/kegg_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$WikiPathway_2023_Human, file = paste(out_path,"/",event_type,"/WIKI/wiki_",cluster,".csv",sep=""), row.names = FALSE)
    write.csv(response_analisi$Reactome_2022, file = paste(out_path,"/",event_type,"/REACTOME/reactome_",cluster,".csv",sep=""), row.names = FALSE)
  }
}
