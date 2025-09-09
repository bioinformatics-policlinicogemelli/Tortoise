"""Module provides various utility functions.

This module provides functions for processing and analyzing
mutational and clinical data, performing enrichment analysis, and creating
and manipulating graphs for visualization and further analysis.

Functions:
    - load_df(config): Load dataframes from CSV files based on the provided
      configuration.
    - filter_vaf(config, df_mut): Filter a DataFrame of mutations based on
      Variant Allele Frequency (VAF) score.
    - enrichment_with_r(path_save, map_cluster): Perform enrichment analysis
      using an R script for given gene clusters.
    - adding_category_mutation(data_mutational, list_columns): Add a category
      mutation label to the mutational data.
    - calculated_vaf(riga): Calculate the Variant Allele Frequency (VAF) if
      absent.
    - create_maps(data_mutational, column_mutation, column_gene, column_sample,
      vaf_score, column_vaf): Create maps for variants and patients.
    - graph_creation(map_patients, map_variants): Create a graph from the maps
      of patients and variants.
    - count_gene(graph): Count the genes present in relation to individual
      mutations.
    - process_data(args): Process data for modularity calculation using the
      Leiden algorithm.
    - selected_seed(g): Select the seed that gives the highest modularity value
      after the Leiden algorithm.
    - leiden_clustering(graph, best_seed): Perform Leiden clustering on the
      graph using the selected seed.
    - adding_graph_color(graph, dendro): Add colors to the graph based on
      clusters for visualization.
    - save_graph_to_file(graph, path_save): Save the graph to a file in GraphML
      format for Cytoscape.
    - map_cluster_creation(graph, dendro): Create a map for clusters from the
      graph and dendrogram.
    - adding_cluster_to_map(map_cluster, map_patients, map_variants): Add
      cluster information to the maps of patients and variants.
    - centroids_cluster(dendro, path_save): Write the centroids of each cluster
      to a file.
    - degree_variant_cluster(map_cluster, graph, path_save): Write the degree
      of each variant in the cluster to a file.
    - enriched_sample_data(data_sample, map_patient, sample_name): Add clinical
      sample information to the patient map.
    - enriched_patient_data(data_patient, map_patient, patient_name): Add
      clinical patient information to the patient map.
    - adding_clinical_info_graph(graph, map_patients): Add clinical information
      to the graph vertices.
    - creation_cluster_clinical_data(map_patient, path_saved): Create a file
      with cluster and clinical data information.
    - numerosity_info(g, map_cluster, path_save): Write a summary of the number
      of patients, variants, and genes in each cluster.
    - summary_info(g, map_cluster, patient_column, path_save): Write a summary
      of mutations and patients for each cluster.
    - count_gene_abs_percent(g, map_cluster, gene_total_count, path_save):
      Count the absolute and percentage presence of genes in each cluster.
    - genes_single_cluster(g, map_cluster, path_save): Create a file for each
      cluster containing the genes present.
    - genes_count_mutation_single_cluster(map_cluster_gene_abs, path_save):
      Create a file for each cluster counting the number of mutations per gene.
    - cluster_noded_attributes(g, map_patient, map_variant): Add cluster
      attributes to the nodes of the graph.
"""

import colorsys
import random
import sys
from multiprocessing import Pool
from pathlib import Path

import igraph as ig
import numpy as np
import pandas as pd
from rpy2 import robjects


def load_df(config):
    """Load dataframes from CSV files based on the provided configuration.

    Args:
        config (dict): A dictionary containing file paths and parameters.
            The dictionary should have the following structure:
            {
                "paths": {
                    # Path to the mutational data CSV file
                    "data_mutational": str,
                    # Separator for the mutational data CSV file
                    "data_mutational_sep": str,
                    # Number of rows to skip in the mutational data CSV file
                    "data_mutational_skip": int,
                    # Path to the clinical sample data CSV file (optional)
                    "data_clinical_sample": str,
                    # Separator for the clinical sample data CSV file
                    "data_clinical_sample_sep": str,
                    # Number of rows to skip in the clinical sample data CSV file
                    "data_clinical_sample_skip": int,
                    # Path to the clinical patient data CSV file (optional)
                    "data_clinical_patient": str,
                    # Separator for the clinical patient data CSV file
                    "data_clinical_patient_sep": str,
                    # Number of rows to skip in the clinical patient data CSV file
                    "data_clinical_patient_skip": int,
                }
            }

    Returns:
        tuple: A tuple containing three pandas DataFrames:
            - df_mut (DataFrame): The mutational data.
            - data_clinical_sample (DataFrame or None): The clinical sample data,
            or None if not provided.
            - data_clinical_patient (DataFrame or None): The clinical patient data,
            or None if not provided.

    """
    df_mut = pd.read_csv(
        config["paths"]["data_mutational"],
        sep=config["paths"]["data_mutational_sep"],
        skiprows=config["paths"]["data_mutational_skip"],
        low_memory=False,
    )
    data_clinical_sample = None
    if config["paths"]["data_clinical_sample"] != "":
        data_clinical_sample = pd.read_csv(
            config["paths"]["data_clinical_sample"],
            sep=config["paths"]["data_clinical_sample_sep"],
            skiprows=config["paths"]["data_clinical_sample_skip"],
            low_memory=False,
        )
    data_clinical_patient = None
    if config["paths"]["data_clinical_patient"] != "":
        data_clinical_patient = pd.read_csv(
            config["paths"]["data_clinical_patient"],
            sep=config["paths"]["data_clinical_patient_sep"],
            skiprows=config["paths"]["data_clinical_patient_skip"],
            low_memory=False,
        )
    return df_mut, data_clinical_sample, data_clinical_patient


def filter_vaf(config, df_mut):
    """Filter a DataFrame of mutations based on Variant Allele Frequency (VAF) score.

    This function filters the mutations in the DataFrame `df_mut` based on the VAF score
    specified in the `config` dictionary. If the `column_vaf` is not specified,
    it defaults to "t_AF" and calculates the VAF.

    Args:
        config (dict): Configuration dictionary containing the following keys:
            - "mutation": A dictionary with the keys:
                - "vaf_score": The VAF score threshold.
                - "vaf_column": The column name for the VAF score in the DataFrame.
            - "paths": A dictionary with the keys:
                - "name": The name of the study.
                - "data_mutational_sep": The separator to use when saving data.
        df_mut (DataFrame): The DataFrame containing mutational data.

    Returns:
        DataFrame: The filtered DataFrame with mutations that meet VAF score threshold.
    Saves:
        A CSV file of the filtered mutations to the path specified in the config file.

    """
    column_vaf = config["mutation"]["vaf_column"]

    if column_vaf != "":
        df_mut = df_mut[df_mut[column_vaf] >= config["mutation"]["vaf_score"]]
    else:
        config["mutation"]["vaf_column"] = "t_AF"
        df_mut["t_AF"] = df_mut.apply(calculated_vaf, axis=1)
        df_mut = df_mut[df_mut["t_AF"] >= config["mutation"]["vaf_score"]]

    df_mut.to_csv(
        Path(
            "study",
            config["name"],
            "input",
            "data_mutational_filtered.txt",
        ),
        index=False,
        sep=config["paths"]["data_mutational_sep"],
    )
    return config, df_mut


def enrichment_with_r(path_save, map_cluster) -> None:
    """Perform enrichment analysis using an R script for given gene clusters.

    This function loads an R script and uses it to perform enrichment analysis
    on gene clusters. The results are saved in specified directories.

    Args:
        path_save (str): The path where the results will be saved.
        map_cluster (dict): A dictionary where keys are cluster identifiers
            and values are lists of genes.

    Returns:
        None

    """
    # Load R script
    robjects.r.source("./lib/enrichment.r")
    # R functions
    r_func = robjects.globalenv["all_analisi"]
    g_path = Path(path_save, "gene_cluster_list")
    o_path = Path(path_save, "pathway_analysis")
    Path(o_path, "GO").mkdir(parents=True, exist_ok=True)
    Path(o_path, "KEGG").mkdir(parents=True, exist_ok=True)
    Path(o_path, "WIKI").mkdir(parents=True, exist_ok=True)
    Path(o_path, "REACTOME").mkdir(parents=True, exist_ok=True)
    if sys.platform.startswith("win") or sys.platform.startswith("linux"):
        with Pool() as p:
            p.map(
                r_func,
                [
                    [
                        str(g_path.joinpath(f"genes_cluster_{c}.csv")),
                        c,
                        str(o_path),
                    ]
                    for c in map_cluster
                ],
            )
    else:
        for c in map_cluster:
            r_func(
                [
                    str(g_path.joinpath(f"genes_cluster_{c}.csv")),
                    c,
                    str(o_path),
                ],
            )


def adding_category_mutation(data_mutational, list_columns):
    nuovi_nomi = {}
    # FILL NA
    for col in list_columns:
        data_mutational = data_mutational.fillna({f"{col}": "N/D"})
    # GROUP
    gruppi_mutazioni = data_mutational.groupby(list_columns)
    # GENERATE ID
    for grp, _ in gruppi_mutazioni:
        # GENE_CHROMOSOME_START_END
        nuovi_nomi[grp] = "Mut_" + "_".join(str(g) for g in grp)
    # ADD COLUMN
    data_mutational["TN_mutation_label"] = data_mutational.apply(
        lambda row: nuovi_nomi[tuple(row[c] for c in list_columns)],
        axis=1,
    )
    return data_mutational


# funzione per cacolare, se assente, la colonna della VAF
def calculated_vaf(riga):
    if riga["t_alt_count"] is np.nan or riga["t_ref_count"] is np.nan:
        return np.nan
    if riga["t_alt_count"] + riga["t_ref_count"] == 0:
        return 0
    return riga["t_alt_count"] / (riga["t_alt_count"] + riga["t_ref_count"])


def create_maps(
    data_mutational,
    column_mutation,
    column_gene,
    column_sample,
    vaf_score,
    column_vaf,
):
    map_variants = {}
    map_patients = {}
    for _, _row in data_mutational.iterrows():
        # GET INFOS FROM DATA
        _paz = str(_row[column_sample])
        _gene = str(_row[column_gene])
        _category = str(_row[column_mutation])

        # se questa colonna è vuota tornerà una stringa vuota
        if "HGVSp_Short" in data_mutational.columns:
            _sost_amm = str(_row["HGVSp_Short"])
        else:
            _sost_amm = ""
        _vaf = float(_row[column_vaf]) if vaf_score else ""

        # creazione dizionario delle varianti
        if _category not in map_variants:
            map_variants[_category] = {}
            map_variants[_category]["patients"] = set()
            map_variants[_category]["sost_amm"] = _sost_amm
            map_variants[_category]["gene"] = _gene
            map_variants[_category]["vaf"] = _vaf

        map_variants[_category]["patients"].add(_paz)

        # creazione dizionario dei pazienti
        if _paz not in map_patients:
            map_patients[_paz] = {}
            map_patients[_paz]["variants"] = set()
        map_patients[_paz]["variants"].add(_category)

    return map_patients, map_variants


# Creazione del Grafo
def graph_creation(map_patients, map_variants):
    edges = [
        (_k_variant, _k_patient)
        for _k_variant, _v_variant in map_variants.items()
        for _k_patient in map_patients
        if _k_patient in _v_variant["patients"]
    ]

    graph = ig.Graph()
    graph.add_vertices(
        list(map_variants.keys()),
        attributes={
            "vertex_type": "VARIANT",
            "color_vertex": "blue",
            "shape_vertex": "circle",
            "gene": [f"{value['gene']}" for value in map_variants.values()],
            "sost_amm": [
                f"{value['sost_amm']}" for value in map_variants.values()
            ],
        },
    )
    graph.add_vertices(
        list(map_patients.keys()),
        attributes={
            "vertex_type": "PATIENT",
            "color_vertex": "red",
            "shape_vertex": "triangle-up",
        },
    )
    graph.add_edges(edges)
    return graph


# count dei geni presenti in relazione alle singole mutazioni
def count_gene(graph):
    gene_total_count = {}
    for vertex in graph.vs:
        if vertex["vertex_type"] == "VARIANT":
            if vertex["gene"] not in gene_total_count:
                gene_total_count[vertex["gene"]] = 0
            gene_total_count[vertex["gene"]] += 1
    return dict(
        sorted(gene_total_count.items(), key=lambda kv: kv[1], reverse=True),
    )


def process_data(args):
    g = args[0]
    _seed = args[1]
    random.seed(_seed)
    np.random.seed(_seed)
    _dendro_2 = g.community_leiden(objective_function="modularity")
    return _dendro_2.modularity


# SELEZIONE DEL SEED CHE Dà VALORE DI MODULARITà PIù ALTA A SEGUITO DEL LEIDEN ALGORITHM
def selected_seed(g, seed_trials):
    mod_results = []
    data = [(g, s) for s in range(seed_trials)]
    if sys.platform.startswith("win") or sys.platform.startswith("linux"):
        with Pool() as p:
            mod_results = p.map(process_data, data)
    else:
        for s in data:
            mod_results.append(process_data(s))

    return mod_results.index(max(mod_results))


# SELEZIONE DEL SEED E LANCIO DELL'ALGORITMO DI CLUSTERIZZAZIONE
def leiden_clustering(graph, best_seed):
    random.seed(best_seed)
    return graph.community_leiden(objective_function="modularity")


# adding color for cluster (for cytoscape)
def adding_graph_color(graph, dendro):
    num_clusters = len(dendro)
    # Genera colori in base al numero di cluster
    colors = []
    for i in range(num_clusters):
        hue = (
            i / num_clusters
        )  # variare la tonalità in base al numero di cluster
        rgb = colorsys.hsv_to_rgb(
            hue,
            1,
            1,
        )  # converti da spazio colore HSV a RGB
        colors.append(
            "#{:02x}{:02x}{:02x}".format(*tuple(int(c * 255) for c in rgb)),
        )  # formato esadecimale RGB
    # Assegna i colori alle singole comunità nel grafo
    graph.vs["color"] = [colors[cluster] for cluster in dendro.membership]

    return graph


# function to save graph as graphml file for cytoscape
def save_graph_to_file(graph, path_save) -> None:
    graph.write_graphml(f"{path_save}/graph_cytoscape.graphml")
    np.save(Path(path_save, "graph.npy"), graph)


# function to create a map for cluster
def map_cluster_creation(graph, dendro):
    map_cluster = {}
    for cluster_index, cluster in enumerate(dendro):
        map_cluster[cluster_index] = []
        for element in cluster:
            _vertex = graph.vs()[element]
            map_cluster[cluster_index].append(
                (_vertex["name"], _vertex["vertex_type"]),
            )
    return map_cluster


# function to add the cluster to patients and variants's map
def adding_cluster_to_map(map_cluster, map_patients, map_variants):
    for cluster, infos in map_cluster.items():
        for info in infos:
            if info[1] == "PATIENT":
                map_patients[info[0]]["cluster"] = cluster
            elif info[1] == "VARIANT":
                map_variants[info[0]]["cluster"] = cluster

    return map_patients, map_variants


# function to write the centroids of each cluster into a file
def centroids_cluster(dendro, path_save) -> None:
    with Path(path_save, "mutation_centroids.csv").open(
        "w",
        encoding="utf-8",
    ) as f:
        for _i, _ in enumerate(dendro):
            sub_graph = dendro.subgraph(_i)
            max_value = 0
            cen_list = []
            for _v in sub_graph.vs():
                if _v["vertex_type"] != "VARIANT":
                    continue
                temp_val = sub_graph.neighborhood_size(_v)
                if temp_val > max_value:
                    max_value = temp_val
                    cen_list = [(_v["name"], max_value - 1)]
                elif temp_val == max_value:
                    cen_list.append((_v["name"], max_value - 1))
            f.write(
                f"Cluster {_i} centroids found {len(cen_list)}: {cen_list}\n",
            )


# funzione per scrivere il numero di connessioni che ciascuna variante ha nel cluster
def degree_variant_cluster(map_cluster, graph, path_save) -> None:
    Path(path_save, "variants_degree").mkdir(parents=True, exist_ok=True)
    for cluster_index in map_cluster:
        list_vertices_filtered = graph.vs.select(
            lambda x, ci=cluster_index: x["cluster"] == ci,
        )
        g_cluster = graph.induced_subgraph(list_vertices_filtered)
        degrees = g_cluster.degree()
        with Path(
            path_save,
            "variants_degree",
            f"variants_degree_cluster_{cluster_index}.csv",
        ).open(
            "w",
            encoding="utf-8",
        ) as f:
            f.write("Variants\tDegree\n")
            for i, degree in enumerate(degrees):
                if g_cluster.vs[i]["vertex_type"] != "PATIENT":
                    f.write(f"{g_cluster.vs[i]['name']}\t{degree}\n")


# AGGIUNTA DELLE INFORMAZIONI CLINICHE DI INTERESSE ALLA MAPPA DEI PAZIENTI
# + AGGIUNTA DEL CLUSTER DI APPARTENENZA ALLA MAPPA DEI PAZIENTI E DELLE MITAZIONI
def enriched_sample_data(data_sample, map_patient, sample_name):
    # generalizzazione --> CLINICAL SAMPLE
    for _, _row in data_sample.iterrows():
        _sample = _row[sample_name]
        for parameter in data_sample.columns:
            variable = _row[parameter]
            if _sample in map_patient:
                map_patient[_sample][parameter] = variable
    return map_patient


def enriched_patient_data(data_patient, map_patient, patient_name):
    # generalizzazione --> CLINICAL PATIENT
    for _, _row in data_patient.iterrows():
        _paz = _row[patient_name]
        for parameters in data_patient.columns:
            variable = _row[parameters]
            for _sample, value in map_patient.items():
                if value[patient_name] == _paz:
                    map_patient[_sample][parameters] = variable
    return map_patient


# function to add clinical information to the vertex
def adding_clinical_info_graph(graph, map_patients):
    for vertex in graph.vs:
        nome_paziente = vertex["name"]
        if nome_paziente in map_patients:
            # Aggiungi le info cliniche come attributi del nodo
            for key, value in map_patients[nome_paziente].items():
                vertex[key] = value
    return graph


# CREAZIONE DEL FILE "CLUSTER_CLINICAL_DATA" IN CUI INSERIAMO
# IL CLUSTER DI APPARTENENZA DI OGI SAMPLES + INFORMAZIONI CLINICHE
def creation_cluster_clinical_data(map_patient, path_saved) -> None:
    header_clinical_data = sorted(
        {k2 for v in map_patient.values() for k2 in v if k2 != "variants"},
    )
    with Path(path_saved, "cluster_clinical_data.csv").open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write("\t".join(header_clinical_data) + "\n")
        for v in map_patient.values():
            temp = [str(v[data]) for data in header_clinical_data]
            f.write("\t".join(temp) + "\n")


# funzione che scrive un file che riassume le info presenti in ogni cluster
def numerosity_info(g, map_cluster, path_save) -> None:
    with Path(path_save, "numerosity_cluster.csv").open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write("Cluster\tPatient\tVariant\tGene\n")
        for cluster in map_cluster:
            count_variant = len(
                [
                    v
                    for v in g.vs
                    if v["cluster"] == cluster
                    and v["vertex_type"] == "VARIANT"
                ],
            )
            count_patient = len(
                [
                    v["gene"]
                    for v in g.vs
                    if v["cluster"] == cluster
                    and v["vertex_type"] == "PATIENT"
                ],
            )
            gene_count = len(
                {
                    v["gene"]
                    for v in g.vs
                    if v["cluster"] == cluster
                    and v["vertex_type"] == "VARIANT"
                },
            )
            f.write(
                str(cluster)
                + "\t"
                + str(count_patient)
                + "\t"
                + str(count_variant)
                + "\t"
                + str(gene_count)
                + "\n",
            )


# riassunto delle informazioni (mutazioni e pazienti per ciascun cluster)
def summary_info(g, map_cluster, patient_column, path_save) -> None:
    with Path(path_save, "summury_file.csv").open("w", encoding="utf-8") as f:
        f.write("Cluster\tNode_Name\tType\tGene\tPatient_id\n")
        for cluster in map_cluster:
            for v in [v for v in g.vs if v["cluster"] == cluster]:
                if v["vertex_type"] == "VARIANT":
                    f.write(
                        str(cluster)
                        + "\t"
                        + str(v["name"])
                        + "\tVARIANT\t"
                        + str(v["gene"])
                        + "\t"
                        + "None"
                        + "\n",
                    )
                else:
                    f.write(
                        str(cluster)
                        + "\t"
                        + str(v["name"])
                        + "\tPATIENT\t"
                        + "None"
                        + "\t",
                    )
                    if patient_column:
                        f.write(str(v[patient_column]))
                    f.write("\n")


# COUNT DEI GENI PRESENTI IN OGNI CLUSTER
# + VALORE PERCENTUALE DI APPARTENENZA DI OGNI GENE A OGNI CLUSTER,
# CON SALVATAGGIO DELLE INFO IN UN FILE "DISTRIBUTION_GENE_CLUSTER"
def count_gene_abs_percent(g, map_cluster, gene_total_count, path_save):
    map_cluster_gene = {}
    for cluster in map_cluster:
        map_cluster_gene[cluster] = {}
        _list_genes = [
            v["gene"]
            for v in g.vs
            if v["cluster"] == cluster and v["vertex_type"] == "VARIANT"
        ]
        for _g in set(_list_genes):
            count = _list_genes.count(_g)
            map_cluster_gene[cluster][_g] = count

    map_cluster_percent = {}
    for cluster, genes in map_cluster_gene.items():
        map_cluster_percent[cluster] = {}
        for gene in genes:
            if gene not in map_cluster_percent[cluster]:
                map_cluster_percent[cluster][gene] = (
                    map_cluster_gene[cluster][gene] * 100
                ) / gene_total_count[gene]

    map_cluster_percent = {
        k: dict(sorted(v.items(), key=lambda x: x[1], reverse=True))
        for k, v in map_cluster_percent.items()
    }
    with Path(path_save, "distribution_gene_cluster.csv").open(
        "w",
        encoding="utf-8",
    ) as f:
        f.write("Cluster\tGene\tPercentage\n")
        for cluster, genes in map_cluster_percent.items():
            for gene, count in genes.items():
                f.write(
                    str(cluster) + "\t" + str(gene) + "\t" + str(count) + "\n",
                )

    return map_cluster_gene, map_cluster_percent


# creazione di un file per ogni cluster, contenente i geni presenti
def genes_single_cluster(g, map_cluster, path_save) -> None:
    Path(path_save, "gene_cluster_list").mkdir(parents=True, exist_ok=True)
    for cluster in map_cluster:
        with Path(path_save, "gene_cluster_list", f"genes_cluster_{cluster}.csv").open(
            "w",
            encoding="utf-8",
        ) as f:
            set_gene = {
                v["gene"]
                for v in g.vs
                if v["cluster"] == cluster and v["vertex_type"] == "VARIANT"
            }
            for genes in set_gene:
                f.write(genes + "\n")


# creazione di un file che per ogni cluster,
# tiene il conto del numero di mutazioni presenti sul gene
def genes_count_mutation_single_cluster(
    map_cluster_gene_abs,
    path_save,
) -> None:
    Path(path_save, "count_cluster_list").mkdir(parents=True, exist_ok=True)
    for cluster, infos in map_cluster_gene_abs.items():
        with Path(
            path_save,
            "count_cluster_list",
            f"count_cluster_{cluster}.csv",
        ).open(
            "w",
            encoding="utf-8",
        ) as f:
            f.write("GENE\tCOUNT\n")
            for k, v in infos.items():
                f.write(f"{k}\t{v}\n")


# AGGIUNTA DELL'ATTRIBUTO CLUSTER AI NODI DEL GRAFO
def cluster_noded_attributes(g, map_patient, map_variant):
    for v in g.vs():
        if v["vertex_type"] == "PATIENT":
            v["cluster"] = map_patient[v["name"]]["cluster"]
        else:
            v["cluster"] = map_variant[v["name"]]["cluster"]
    return g
