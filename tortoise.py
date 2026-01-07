#!/usr/bin/env python3
"""Module that creates graphs and performs clustering and enrichment analysis.

This module processes mutational and clinical data, creates bipartite graphs
(patients-variants for mutations and patients-genes for CNV), performs Leiden
clustering, and analyzes gene/CNV enrichment by cluster.

Functions:
    main(path_config): Main function to process data and perform analysis.
        Handles both mutation and optional CNV data to create and analyze graphs.
    load_df(config): Loads mutational and clinical data from CSV files.
    filter_vaf(config, df_mut): Filters mutations based on VAF score.
    enrichment_with_r(path_save, map_c): Performs enrichment analysis.

Features:
    - Mutation graph analysis: Creates bipartite graphs of patients and variants
    - CNV analysis: Optional Copy Number Variation analysis with separate graph
    - Leiden clustering: Community detection with seed optimization
    - Enrichment analysis: Gene set enrichment using R scripts
    - Output generation: CSV files with cluster statistics and centroids

Usage:
    Run this module as a script with the required configuration file path:
    python tortoise.py -c /path/to/config.json
"""

import argparse
import json
import logging
from pathlib import Path

import lib.lib_utils as libu

logger = logging.getLogger("tortoise")
logging.basicConfig(
    format="%(asctime)s %(name)s -- %(levelname)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    encoding="utf-8",
)


def main(path_config: Path) -> None:
    """Create graphs, and perform clustering and enrichment analysis.

    Args:
        path_config (Path): Path to the configuration file in JSON format.
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
    logger.info("  0% -- Start tortoise")
    with Path(path_config).open(mode="r", encoding="utf-8") as f:
        config = json.load(f)
    path_save = Path("study", config["name"], "output")
    identifier_columns = config["mutation"]["identifier_columns"].split(";")
    column_mutation_name = None
    Path(path_save).mkdir(parents=True, exist_ok=True)
    # load df
    logger.info(" 10% -- Load files")
    df_mut, data_clinical_sample, data_clinical_patient, data_cnv = libu.load_df(config)
    # create maps and graph
    logger.info(" 20% -- Add category mutations")
    if len(identifier_columns) > 1:
        df_mut = libu.adding_category_mutation(df_mut, identifier_columns)
        column_mutation_name = "TN_mutation_label"
    else:
        column_mutation_name = identifier_columns[0]
    # VAF filter
    if config["mutation"]["vaf_score"]:
        config, df_mut = libu.filter_vaf(config, df_mut)
    logger.info(" 25% -- Create maps")
    map_patients, map_variants = libu.create_maps(
        df_mut,
        column_mutation_name,
        config["mutation"]["column_gene"],
        config["mutation"]["column_sample_name"],
        config["mutation"]["vaf_score"],
        config["mutation"]["vaf_column"],
    )
    # Add CNV
    logger.info(" 30% -- Add CNV analysis if data is available")
    cnv_available = not (data_cnv is None or data_cnv.empty)
    if cnv_available:
        logger.info(" CNV detected: processing...")

        data_cnv = libu.preprocess_cnv(data_cnv)

        cnv_identifier_columns = config["cnv"]["column_cnv_identifier"].split(";")

        if len(cnv_identifier_columns) > 1:
            data_cnv = libu.adding_category_mutation(data_cnv, cnv_identifier_columns)
            column_cnv_name = "TN_cnv_label"
        else:
            column_cnv_name = cnv_identifier_columns[0]

        map_pat, map_cnv = libu.create_cnv_map(data_cnv)
    else:
        logger.info(" No CNV file provided: skipping CNV analysis.")
        data_cnv = None
        map_pat, map_cnv = None, None

    # cluster
    logger.info(" 40% -- Create graph")
    g = libu.graph_creation(map_patients, map_variants)
    logger.info(" 50% -- Clustering")
    best_seed = libu.selected_seed(g, config["seed_trials"], config["clustering_resolution"])
    dendro = libu.leiden_clustering(g, best_seed, config["clustering_resolution"])
    with Path(path_save, "modularity.info").open("w") as f:
        f.write(str(round(dendro.modularity, 4)))
    with Path(path_save, "seed.info").open("w") as f:
        f.write(str(best_seed))
    logger.info(f" best seed: {best_seed} -- modularity: {round(dendro.modularity, 4)}")
    # graph add colors + cytoscape
    g = libu.adding_graph_color(g, dendro)
    # map cluster
    logger.info(" 60% -- Create map cluster")
    map_c = libu.map_cluster_creation(g, dendro)
    map_patients, map_variants = libu.adding_cluster_to_map(
        map_c,
        map_patients,
        map_variants,
    )
    # cluster CNV ############################### IF CNV AVAILABLE
    if cnv_available:
        logger.info(" 63% -- Create CNV graph")
        g2 = libu.graph_cnv(map_pat, map_cnv)

        logger.info(" 66% -- CNV Clustering")
        bs2 = libu.selected_seed(g2, config["seed_trials"], config["clustering_resolution"])
        dendro_cnv = libu.leiden_clustering(g2, bs2, config["clustering_resolution"])

        with Path(path_save, "modularity_cnv.info").open("w") as f:
            f.write(str(round(dendro_cnv.modularity, 4)))

        with Path(path_save, "seed_cnv.info").open("w") as f:
            f.write(str(bs2))

        # colors, mapping, attributes
        g2 = libu.adding_graph_color(g2, dendro_cnv)
        map_c2 = libu.map_cluster_creation(g2, dendro_cnv)
        map_pat, map_cnv = libu.adding_cluster_to_map_cnv(map_c2, map_pat, map_cnv)

        g2 = libu.cluster_cnv_noded_attributes(g2, map_pat, map_cnv)
        libu.save_graph_to_file(g2, path_save, name="graph_cnv")
    else:
        g2 = None
        dendro_cnv = None
        map_c2 = None
    logger.info(" 70% -- Add metadata")
    g = libu.cluster_noded_attributes(g, map_patients, map_variants)
    # add info
    if data_clinical_sample is not None:
        map_patients = libu.enriched_sample_data(
            data_clinical_sample,
            map_patients,
            config["clinical_data"]["column_sample_name"],
        )
    if data_clinical_patient is not None:
        map_patients = libu.enriched_patient_data(
            data_clinical_patient,
            map_patients,
            config["clinical_data"]["column_patient_name"],
        )
    g = libu.adding_clinical_info_graph(g, map_patients)
    logger.info(" 80% -- Export infos")
    libu.save_graph_to_file(g, path_save, name="graph_mutational")
    # create info files
    libu.summary_info(
        g,
        map_c,
        config["clinical_data"]["column_patient_name"],
        path_save,
    )
    libu.numerosity_info(g, map_c, path_save)
    # stats
    map_cluster_gene_abs, _ = libu.count_gene_abs_percent(
        g,
        map_c,
        libu.count_gene(g),
        path_save,
    )
    libu.genes_single_cluster(g, map_c, path_save)
    libu.genes_count_mutation_single_cluster(map_cluster_gene_abs, path_save)
    libu.creation_cluster_clinical_data(map_patients, path_save)
    libu.centroids_cluster(dendro, path_save)
    libu.degree_variant_cluster(map_c, g, path_save)

    if cnv_available:
        libu.numerosity_info_cnv(g2, map_c2, path_save)
        map_cluster_cnv_abs, _ = libu.count_gene_abs_percent_cnv(
            g2,
            map_c2,
            libu.count_gene_cnv(g2),
            path_save,
        )
        libu.genes_single_cluster_cnv(g2, map_c2, path_save)
        libu.genes_count_cnv_single_cluster(map_cluster_cnv_abs, path_save)
        libu.degree_variant_cluster_cnv(map_c2, g2, path_save)
        libu.centroids_cluster_cnv(dendro_cnv, path_save)

    # enrichment
    logger.info(" 90% -- Enrichment")
    libu.enrichment_with_r(path_save, map_c)
    logger.info("100% -- Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--verbose",
        help="increase verbosity",
        action="store_true",
    )
    parser.add_argument("-c", "--config", type=str, required=True)
    main(Path(parser.parse_args().config))
