#!/usr/bin/env python3
"""Module taht creates graphs and performs clustering and enrichment analysis.

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
    df_mut, data_clinical_sample, data_clinical_patient = libu.load_df(config)
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
    logger.info(" 30% -- Create maps")
    map_patients, map_variants = libu.create_maps(
        df_mut,
        column_mutation_name,
        config["mutation"]["column_gene"],
        config["mutation"]["column_sample_name"],
        config["mutation"]["vaf_score"],
        config["mutation"]["vaf_column"],
    )
    # cluster
    logger.info(" 40% -- Create graph")
    g = libu.graph_creation(map_patients, map_variants)
    logger.info(" 50% -- Clustering")
    dendro = libu.leiden_clustering(g, libu.selected_seed(g))
    with Path(path_save, "modularity.info").open("w") as f:
        f.write(str(round(dendro.modularity, 4)))
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
    libu.save_graph_to_file(g, path_save)
    # create inofo files
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
