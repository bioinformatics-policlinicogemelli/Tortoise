import os
import sys
import json
import argparse
import pandas as pd
import lib.lib_utils as libu
import rpy2.robjects as robjects
from multiprocessing import Pool

def exec(path_config):
    #CARICAMENTO FILE CONFIGURAZIONE
    with open(path_config) as f:
        DATA = json.load(f)
    #Definizione pathway
    path_save                     = os.path.join("study",DATA["Paths"]["name_study"],"output")
    column_mutation_name          = DATA["Mutation"]["column_mutation_name"]
    column_gene                   = DATA["Mutation"]["column_gene_name"]
    column_hgvsp                  = DATA["Mutation"]["column_hgvsp_short"]
    column_variant_classification = DATA["Mutation"]["column_variant_classification"]
    column_chromosome             = DATA["Mutation"]["column_chromosome"]
    column_start                  = DATA["Mutation"]["column_start"]
    column_end                    = DATA["Mutation"]["column_end"]
    sample_name                   = DATA["Clinical_data"]["column_sample_name"]
    patient_name                  = DATA["Clinical_data"]["column_patient_name"]
    vaf                           = DATA["Mutation"]["vaf"]
    vaf_score                     = DATA["Mutation"]["vaf_score"]
    column_vaf                    = DATA["Mutation"]["vaf_column"]
    os.makedirs(path_save, exist_ok=True)
    #caricamento dei dataframes
    data_mutational = pd.read_csv(
        DATA["Paths"]["data_mutational"],
        sep = DATA["Paths"]["data_mutational_sep"],
        skiprows = DATA["Paths"]["data_mutational_skip"],
        low_memory = False
    )
    data_clinical_sample = None
    if DATA["Paths"]["data_clinical_sample"] != "":
        data_clinical_sample = pd.read_csv(
        DATA["Paths"]["data_clinical_sample"],
        sep = DATA["Paths"]["data_clinical_sample_sep"],
        skiprows = DATA["Paths"]["data_clinical_sample_skip"],
        low_memory = False
    )
    data_clinical_patient = None
    if DATA["Paths"]["data_clinical_patient"] != "":
        data_clinical_patient = pd.read_csv(
        DATA["Paths"]["data_clinical_patient"],
        sep = DATA["Paths"]["data_clinical_patient_sep"],
        skiprows = DATA["Paths"]["data_clinical_patient_skip"],
        low_memory = False
    )
    # creazione delle mappe pazienti e varianti + creazione del grafo
    if column_mutation_name == "":
        data_mutational = libu.adding_category_mutation(
            data_mutational,
            column_gene,
            column_hgvsp,
            column_variant_classification,
            column_chromosome,
            column_start,
            column_end
        )
    #filtraggio della VAF
    if vaf:
        if column_vaf != "":
            data_mutational = data_mutational[data_mutational[column_vaf] >= vaf_score]
        else:
            data_mutational['t_AF'] = data_mutational.apply(libu.calculated_vaf, axis = 1)
            data_mutational = data_mutational[data_mutational['t_AF'] >= vaf_score]

        data_mutational.to_csv(
            os.path.join("study",DATA["Paths"]["name_study"],"input","data_mutational_filtered.txt"),
            index = False,
            sep = DATA["Paths"]["data_mutational_sep"]
        )
    map_patients,map_variants = libu.create_maps(data_mutational,DATA,column_mutation_name)
    #clusterizzazione
    graph=libu.graph_creation(map_patients,map_variants)
    seed=libu.selected_seed(graph)
    dendro=libu.leiden_clustering(graph,seed)
    #gestione parte grafica del grafo (aggiunta colori + file per cytoscape)
    graph=libu.adding_graph_color(graph,dendro)
    #creazione della mappa cluster e attribuzione del cluster i pazienti e alle varianti
    map_cluster=libu.map_cluster_creation(graph,dendro)
    map_patients,map_variants=libu.adding_cluster_to_map(map_cluster,map_patients,map_variants)
    graph=libu.cluster_noded_attributes(graph,map_patients,map_variants)
    #aggiunta delle informazioni cliniche alla mappa dei pazienti
    if data_clinical_sample is not None:
        map_patients=libu.enriched_sample_data(data_clinical_sample,map_patients,sample_name)
    if data_clinical_patient is not None:
        map_patients=libu.enriched_patient_data(data_clinical_patient,map_patients,patient_name)
    graph=libu.adding_clinical_info_graph(graph,map_patients)
    libu.save_graph_to_file(graph, path_save)
    #creazione file in cui riassumere le informazioni nei diversi cluster
    libu.summary_info(path_save,map_cluster,map_patients,patient_name)
    libu.numerosity_info(path_save,map_cluster)
    #creazione di una mappa con il numero di mutazioni per ogni gene + creazione di due mappe con i valori assoluti e 
    #percentuali di distribuzione delle mutazioni, per ciascun gene,nei diversi cluster
    gene_total_count=libu.count_gene(graph)
    #FIXME rivedi funzione
    map_cluster_gene_abs,map_cluster_gene_percent=libu.count_gene_abs_percent(map_cluster,gene_total_count,path_save)
    #salvataggio delle percentuali di distribuzione delle mutazioni dei diversi geni nei cluster
    libu.genes_single_cluster(map_cluster,path_save)
    libu.genes_count_mutation_single_cluster(map_cluster_gene_abs,path_save)
    libu.creation_cluster_clinical_data(map_patients,path_save)
    libu.centroids_cluster(dendro,path_save)
    libu.degree_variant_cluster(map_cluster,graph,path_save)
    # Arricchimento
    #Load R script
    robjects.r.source("./lib/enrichment.r")
    #R functions
    all_analisi             = robjects.globalenv['all_analisi']
    genes_path = os.path.join(path_save, "Gene")
    out_path = os.path.join(path_save,"Arricchimento_all_genes")
    os.makedirs(os.path.join(out_path,"GO"), exist_ok=True)
    os.makedirs(os.path.join(out_path,"KEGG"), exist_ok=True)
    os.makedirs(os.path.join(out_path,"WIKI"), exist_ok=True)
    os.makedirs(os.path.join(out_path,"REACTOME"), exist_ok=True)
    if sys.platform.startswith('win') or sys.platform.startswith("linux"):
        with Pool() as p:
            p.map(all_analisi, [[os.path.join(genes_path,f"genes_cluster_{c}.csv"),c,out_path] for c in map_cluster.keys()])
    else:
        for c in map_cluster.keys():
            all_analisi([os.path.join(genes_path,f"genes_cluster_{c}.csv"),c,out_path])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-c", "--config", type=str, required=True)
    args = parser.parse_args()
    exec(args.config)
