"""
This module processes mutational and clinical data,
creates graphs, and performs clustering and enrichment analysis.

Functions:
    load_df(config): Loads mutational and clinical data from CSV files.
    filter_vaf(config, df_mut): Filters mutations based on VAF score.
    enrichment_with_r(path_save, map_c): Performs enrichment analysis.
    main(path_config): Main function to process data and perform analysis.

Usage:
    Run this module as a script with the required configuration file path:
    python tortoise.py -c /path/to/config.json
"""

import argparse
import json
import os

import lib.lib_utils as libu


def main(path_config):
    """
    Main function to process mutational and clinical data,
    create graphs, and perform clustering and enrichment analysis.
    Args:
        path_config (str): Path to the configuration file in JSON format.
    The function performs the following steps:
    1. Loads the configuration file.
    2. Defines paths and loads mutational and clinical data.
    3. Creates patient and variant maps and constructs a g.
    4. Filters the data based on Variant Allele Frequency (VAF).
    5. Performs clustering using the Leiden algorithm.
    6. Adds graphical attributes and clinical information to the g.
    7. Saves the g and generates summary information.
    8. Counts gene mutations and creates distribution maps.
    9. Performs enrichment analysis using an R script.
    """

    # load config file
    with open(path_config, mode="r", encoding="utf-8") as f:
        config = json.load(f)
    path_save = os.path.join("study", config["Paths"]["name_study"], "output")
    identifier_columns = config["Mutation"]["identifier_columns"].split(";")
    column_mutation_name = None
    os.makedirs(path_save, exist_ok=True)
    # load df
    df_mut, data_clinical_sample, data_clinical_patient = libu.load_df(config)
    # creazione delle mappe pazienti e varianti + creazione del grafo
    if len(identifier_columns) > 1:
        df_mut = libu.adding_category_mutation(df_mut, identifier_columns)
        column_mutation_name = "TN_mutation_label"
    else:
        column_mutation_name = identifier_columns[0]
    # filtraggio della VAF
    if config["Mutation"]["vaf"]:
        config, df_mut = libu.filter_vaf(config, df_mut)
    map_patients, map_variants = libu.create_maps(
        df_mut,
        column_mutation_name,
        config["Mutation"]["column_gene"],
        config["Mutation"]["column_sample_name"],
        config["Mutation"]["vaf"],
        config["Mutation"]["vaf_column"],
    )
    # clusterizzazione
    g = libu.graph_creation(map_patients, map_variants)
    dendro = libu.leiden_clustering(g, libu.selected_seed(g))
    # gestione parte grafica del grafo (aggiunta colori + file per cytoscape)
    g = libu.adding_graph_color(g, dendro)
    # creazione della mappa cluster
    map_c = libu.map_cluster_creation(g, dendro)
    map_patients, map_variants = libu.adding_cluster_to_map(
        map_c, map_patients, map_variants
    )
    g = libu.cluster_noded_attributes(g, map_patients, map_variants)
    # aggiunta delle informazioni cliniche alla mappa dei pazienti
    if data_clinical_sample is not None:
        map_patients = libu.enriched_sample_data(
            data_clinical_sample,
            map_patients,
            config["Clinical_data"]["column_sample_name"],
        )
    if data_clinical_patient is not None:
        map_patients = libu.enriched_patient_data(
            data_clinical_patient,
            map_patients,
            config["Clinical_data"]["column_patient_name"],
        )
    g = libu.adding_clinical_info_graph(g, map_patients)
    libu.save_graph_to_file(g, path_save)
    # creazione file in cui riassumere le informazioni
    libu.summary_info(
        g, map_c, config["Clinical_data"]["column_patient_name"], path_save
    )
    libu.numerosity_info(g, map_c, path_save)
    # creazione di una mappa con il numero di mutazioni per ogni gene
    # creazione di due mappe con i valori assoluti e
    # percentuali di distribuzione delle mutazioni, per ciascun gene
    map_cluster_gene_abs, _ = libu.count_gene_abs_percent(
        g, map_c, libu.count_gene(g), path_save
    )
    # percentuali di distribuzione delle mutazioni dei diversi geni nei cluster
    libu.genes_single_cluster(g, map_c, path_save)
    libu.genes_count_mutation_single_cluster(map_cluster_gene_abs, path_save)
    libu.creation_cluster_clinical_data(map_patients, path_save)
    libu.centroids_cluster(dendro, path_save)
    libu.degree_variant_cluster(map_c, g, path_save)
    libu.enrichment_with_r(path_save, map_c)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--verbose", help="increase verbosity", action="store_true"
    )
    parser.add_argument("-c", "--config", type=str, required=True)
    main(parser.parse_args().config)
