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
import re
import igraph as ig
import numpy as np
import pandas as pd
from rpy2 import robjects

# ============================================================
# MULTI-RESOLUTION CLUSTERING (ADD-ON, NON USATO DI DEFAULT)
# ============================================================

def multi_resolution_clustering(
    graph,
    resolutions,
    seed_trials,
    path_save,
):
    """
    Perform Leiden clustering at multiple resolutions and extract
    centroid clusters.

    For each resolution:
    - selects best seed (max modularity)
    - performs Leiden clustering
    - identifies centroid cluster (highest-degree VARIANT)
    - extracts PATIENT sample IDs
    - writes a TSV with one column per resolution

    Output:
        resolution_centroids.tsv
    """

    centroid_samples = {}

    for res in resolutions:
        print(f"[RESOLUTION] Processing resolution {res}")

        # ---------------------------
        # SEED SELECTION
        # ---------------------------
        def _process(args):
            g, seed, resolution = args
            random.seed(seed)
            np.random.seed(seed)
            d = g.community_leiden(
                objective_function="modularity",
                resolution_parameter=resolution,
            )
            return d.modularity

        data = [(graph, s, res) for s in range(seed_trials)]

        if sys.platform.startswith(("linux", "win")):
            with Pool() as p:
                modularities = p.map(_process, data)
        else:
            modularities = [_process(d) for d in data]

        best_seed = modularities.index(max(modularities))

        # ---------------------------
        # FINAL CLUSTERING
        # ---------------------------
        random.seed(best_seed)
        dendro = graph.community_leiden(
            objective_function="modularity",
            resolution_parameter=res,
        )

        # ---------------------------
        # PER-CLUSTER CENTROIDS
        # ---------------------------
        for cluster_id in range(len(dendro)):
            sub = dendro.subgraph(cluster_id)

            # calcolo degree massimo tra le VARIANT
            max_degree = -1
            centroid_genes = set()

            for v in sub.vs:
                if v["vertex_type"] == "VARIANT":
                    deg = sub.degree(v.index)
                    if deg > max_degree:
                        max_degree = deg
                        centroid_genes = {v["gene"]}
                    elif deg == max_degree:
                        centroid_genes.add(v["gene"])

            if not centroid_genes:
                continue  # cluster senza varianti (caso raro)

            # sampleID del cluster
            samples = [
                v["name"]
                for v in sub.vs
                if v["vertex_type"] == "PATIENT"
            ]

            if not samples:
                continue

            # header colonna
            gene_label = "-".join(sorted(centroid_genes))
            column_name = f"{gene_label}_cluster{cluster_id}_{res}"

            centroid_samples[column_name] = sorted(samples)

    # ---------------------------
    # WRITE TSV (ALL RESOLUTIONS)
    # ---------------------------
    max_len = max(len(v) for v in centroid_samples.values())

    df = pd.DataFrame(
        {
            k: v + [""] * (max_len - len(v))
            for k, v in centroid_samples.items()
        }
    )

    out_file = Path(path_save, "resolution_centroids.tsv")
    df.to_csv(out_file, sep="\t", index=False)

    print(f"[DONE] Resolution analysis saved to {out_file}")

# =================================
# NORMALIZATION DATA SURVIVAL EVENT
# =================================

def normalize_survival_event(series: pd.Series) -> pd.Series:
    """
    Normalize survival event column to:
    - 1.0 = event occurred
    - 0.0 = censored
    - NaN = unknown / invalid

    Handles:
    - 1 / 0
    - yes / no
    - dead / alive
    - progress / progressed / progression
    - censored / censor
    - compound values like '1:PROGRESS', '0:CENSORED'
    - case-insensitive, punctuation-agnostic
    """

    if series is None:
        return series

    # keywords indicating event or censoring
    event_values = {
        "1", "yes", "y", "true",
        "dead", "death", "deceased",
        "progress", "progressed", "progression",
        "relapse", "event",
    }

    censored_values = {
        "0", "no", "n", "false",
        "alive", "censored", "censor",
    }

    def _normalize_value(val):

        # ----------------------
        # NULL
        # ----------------------
        if pd.isna(val):
            return np.nan

        # ----------------------
        # REAL NUMBERS
        # ----------------------
        if isinstance(val, (int, float)):

            if np.isnan(val):
                return np.nan

            if val == 1:
                return 1.0

            if val == 0:
                return 0.0

            return np.nan

        # ----------------------
        # STRING
        # ----------------------
        s = str(val).strip().lower()

        if not s:
            return np.nan

        # Try direct numeric conversion first
        try:
            f = float(s)

            if f == 1.0:
                return 1.0

            if f == 0.0:
                return 0.0

        except Exception:
            pass

        # Replace punctuation with spaces
        s = re.sub(r"[^\w]+", " ", s)

        tokens = s.split()

        # Event if ANY token matches
        if any(tok in event_values for tok in tokens):
            return 1.0

        # Censored if ANY token matches
        if any(tok in censored_values for tok in tokens):
            return 0.0

        return np.nan

    # Apply + FORCE FLOAT dtype
    out = series.apply(_normalize_value)

    return out.astype("float64")

# ===========================
# READ AND SET SURVIVAL EVENT
# ===========================

def prepare_survival_columns(df_clinical: pd.DataFrame, survival_cfg: dict):
    """
    Prepare canonical survival columns (OS / PFS) in clinical dataframe.

    Uses ONLY explicit configuration from config.json.
    No column name inference is performed.

    Returns
    -------
    df_clinical : pd.DataFrame
        Updated dataframe with OS_TIME, OS_EVENT, PFS_TIME, PFS_EVENT (if available)
    survival_info : dict
        {
            "os_available": bool,
            "pfs_available": bool
        }
    """

    df = df_clinical.copy()

    survival_info = {
        "os_available": False,
        "pfs_available": False,
    }

    if not isinstance(survival_cfg, dict):
        return df, survival_info

    for endpoint in ("os", "pfs"):
        cfg = survival_cfg.get(endpoint, {})
        if not cfg:
            continue

        time_col = cfg.get("time_column")
        event_col = cfg.get("event_column")

        # Both columns must be explicitly defined and present
        if (
            not time_col
            or not event_col
            or time_col not in df.columns
            or event_col not in df.columns
        ):
            continue

        time_series = pd.to_numeric(df[time_col], errors="coerce")
        event_series = normalize_survival_event(df[event_col])

        valid_mask = time_series.notna() & event_series.notna()
        if valid_mask.sum() < 2:
            continue

        df[f"{endpoint.upper()}_TIME"] = time_series
        df[f"{endpoint.upper()}_EVENT"] = event_series
        survival_info[f"{endpoint}_available"] = True

    return df, survival_info

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
    _graph = args[0]
    _seed = args[1]
    _resolution = args[2]
    random.seed(_seed)
    np.random.seed(_seed)
    _dendro_2 = _graph.community_leiden(objective_function="modularity", resolution_parameter=_resolution)
    return _dendro_2.modularity


# SELEZIONE DEL SEED CHE Dà VALORE DI MODULARITà PIù ALTA A SEGUITO DEL LEIDEN ALGORITHM
def selected_seed(graph, seed_trials, resolution):
    mod_results = []
    data = [(graph, s, resolution) for s in range(seed_trials)]
    if sys.platform.startswith("win") or sys.platform.startswith("linux"):
        with Pool() as p:
            mod_results = p.map(process_data, data)
    else:
        for args in data:
            mod_results.append(process_data(args))

    return mod_results.index(max(mod_results))


# SELEZIONE DEL SEED E LANCIO DELL'ALGORITMO DI CLUSTERIZZAZIONE
def leiden_clustering(graph, seed, resolution):
    random.seed(seed)
    return graph.community_leiden(objective_function="modularity", resolution_parameter=resolution)


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


def gene_centroid_stability(path_save):
    """
    Computes gene centroid stability across resolutions.

    Output:
        gene_centroid_stability.tsv
    """

    path = Path(path_save, "resolution_centroids.tsv")
    df = pd.read_csv(path, sep="\t")

    gene_counts = {}
    resolution_set = set()

    for col in df.columns:
        # esempio: TP53-BRCA1_cluster2_0.6
        gene_part, _, res = col.rpartition("_")
        genes = gene_part.split("_")[0].split("-")
        resolution_set.add(res)

        for g in genes:
            gene_counts[g] = gene_counts.get(g, 0) + 1

    total_resolutions = len(resolution_set)

    out = pd.DataFrame(
        {
            "gene": gene_counts.keys(),
            "stability": [
                gene_counts[g] / total_resolutions for g in gene_counts
            ],
        }
    ).sort_values("stability", ascending=False)

    out.to_csv(
        Path(path_save, "gene_centroid_stability.tsv"),
        sep="\t",
        index=False,
    )


def centroid_overlap_matrix(path_save):
    """
    Computes Jaccard overlap between centroid sample sets
    across resolutions.

    Output:
        centroid_overlap.tsv
    """

    df = pd.read_csv(
        Path(path_save, "resolution_centroids.tsv"),
        sep="\t",
    )

    # mappa: resolution -> set(sampleID)
    res_map = {}

    for col in df.columns:
        res = col.rpartition("_")[2]
        samples = set(df[col].dropna()) - {""}
        res_map.setdefault(res, set()).update(samples)

    resolutions = sorted(res_map.keys())
    matrix = []

    for r1 in resolutions:
        row = []
        for r2 in resolutions:
            inter = len(res_map[r1] & res_map[r2])
            union = len(res_map[r1] | res_map[r2])
            row.append(inter / union if union else 0)
        matrix.append(row)

    out = pd.DataFrame(matrix, index=resolutions, columns=resolutions)
    out.to_csv(
        Path(path_save, "centroid_overlap.tsv"),
        sep="\t",
    )


def select_best_resolution(path_save):
    """
    Selects best resolution based on:
    - centroid stability
    - overlap plateau

    Output:
        best_resolution.info
    """

    stability = pd.read_csv(
        Path(path_save, "gene_centroid_stability.tsv"),
        sep="\t",
    )

    overlap = pd.read_csv(
        Path(path_save, "centroid_overlap.tsv"),
        sep="\t",
        index_col=0,
    )

    # score = media overlap + media stability
    overlap_score = overlap.mean(axis=1)
    stability_score = stability.groupby("gene")["stability"].mean().mean()

    scores = overlap_score + stability_score

    best_res = scores.idxmax()

    with open(Path(path_save, "best_resolution.info"), "w") as f:
        f.write(str(best_res))

def gene_centroid_resolution_heatmap(path_save):
    """
    Creates a gene-centroid vs resolution heatmap.

    Output:
        gene_centroid_resolution.tsv
        gene_centroid_resolution_heatmap.png
    """
    import logging
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    import matplotlib.pyplot as plt
    import seaborn as sns

    # -------------------------------------------------
    # Load centroid table
    # -------------------------------------------------
    df = pd.read_csv(
        Path(path_save, "resolution_centroids.tsv"),
        sep="\t",
    )

    # -------------------------------------------------
    # Parse gene and resolution from column names
    # -------------------------------------------------
    records = []

    for col in df.columns:
        # es: TP53-BRCA1_cluster2_0.6
        left, _, res = col.rpartition("_")
        gene_part = left.split("_cluster")[0]
        genes = gene_part.split("-")

        for g in genes:
            records.append(
                {
                    "gene": g,
                    "resolution": float(res),
                    "count": 1,
                }
            )

    long_df = pd.DataFrame(records)

    # -------------------------------------------------
    # Aggregate and normalize
    # -------------------------------------------------
    heatmap_df = (
        long_df
        .groupby(["gene", "resolution"])
        .count()
        .reset_index()
        .pivot(index="gene", columns="resolution", values="count")
        .fillna(0)
    )

    # normalize per gene
    heatmap_norm = heatmap_df.div(
        heatmap_df.max(axis=1),
        axis=0,
    )

    # -------------------------------------------------
    # Save TSV
    # -------------------------------------------------
    heatmap_norm.to_csv(
        Path(path_save, "gene_centroid_resolution.tsv"),
        sep="\t",
    )

    # -------------------------------------------------
    # Plot heatmap
    # -------------------------------------------------
    plt.figure(figsize=(0.6 * heatmap_norm.shape[1] + 4,
                        0.3 * heatmap_norm.shape[0] + 4))

    sns.heatmap(
        heatmap_norm,
        cmap="viridis",
        linewidths=0.5,
        linecolor="grey",
        cbar_kws={"label": "Normalized centroid frequency"},
    )

    plt.xlabel("Resolution")
    plt.ylabel("Gene centroid")
    plt.title("Gene-centroid stability across resolutions")

    plt.tight_layout()
    plt.savefig(
        Path(path_save, "gene_centroid_resolution_heatmap.png"),
        dpi=300,
    )
    plt.close()


def cluster_centroid_resolution_heatmap(path_save):
    """
    Computes cluster-level centroid stability across resolutions
    using Jaccard overlap with previous resolution.

    Output:
        cluster_centroid_resolution.tsv
        cluster_centroid_resolution_heatmap.png
    """

    import matplotlib.pyplot as plt
    import seaborn as sns

    df = pd.read_csv(
        Path(path_save, "resolution_centroids.tsv"),
        sep="\t",
    )

    # ---------------------------------------
    # Parse cluster -> resolution -> samples
    # ---------------------------------------
    cluster_map = {}

    for col in df.columns:
        # gene_clusterX_res
        left, _, res = col.rpartition("_")
        cluster_id = left.split("_cluster")[-1]
        key = f"cluster{cluster_id}"

        samples = set(df[col].dropna()) - {""}

        cluster_map.setdefault(key, {})[float(res)] = samples

    resolutions = sorted(
        {float(c.rpartition("_")[2]) for c in df.columns}
    )

    # ---------------------------------------
    # Jaccard overlap vs previous resolution
    # ---------------------------------------
    records = []

    for cluster, res_map in cluster_map.items():
        for i in range(1, len(resolutions)):
            r1 = resolutions[i - 1]
            r2 = resolutions[i]

            if r1 in res_map and r2 in res_map:
                s1 = res_map[r1]
                s2 = res_map[r2]
                j = len(s1 & s2) / len(s1 | s2) if s1 | s2 else 0
            else:
                j = 0

            records.append(
                {
                    "cluster": cluster,
                    "resolution": r2,
                    "jaccard": j,
                }
            )

    heat_df = (
        pd.DataFrame(records)
        .pivot(index="cluster", columns="resolution", values="jaccard")
        .fillna(0)
    )

    heat_df.to_csv(
        Path(path_save, "cluster_centroid_resolution.tsv"),
        sep="\t",
    )

    # ---------------------------------------
    # Plot
    # ---------------------------------------
    plt.figure(figsize=(0.6 * heat_df.shape[1] + 4,
                        0.3 * heat_df.shape[0] + 4))

    sns.heatmap(
        heat_df,
        cmap="magma",
        linewidths=0.3,
        cbar_kws={"label": "Jaccard overlap"},
    )

    plt.xlabel("Resolution")
    plt.ylabel("Cluster centroid")
    plt.title("Cluster-level centroid stability")

    plt.tight_layout()
    plt.savefig(
        Path(path_save, "cluster_centroid_resolution_heatmap.png"),
        dpi=300,
    )
    plt.close()


def resolution_vs_centroid_gene_count(path_save):
    """
    Computes number of distinct centroid genes per resolution.

    Output:
        resolution_vs_centroid_genes.tsv
        resolution_vs_centroid_genes.png
    """

    import matplotlib.pyplot as plt

    df = pd.read_csv(
        Path(path_save, "resolution_centroids.tsv"),
        sep="\t",
    )

    res_gene_map = {}

    for col in df.columns:
        left, _, res = col.rpartition("_")
        genes = left.split("_cluster")[0].split("-")

        res_gene_map.setdefault(float(res), set()).update(genes)

    out = pd.DataFrame(
        {
            "resolution": sorted(res_gene_map),
            "n_centroid_genes": [
                len(res_gene_map[r]) for r in sorted(res_gene_map)
            ],
        }
    )

    out.to_csv(
        Path(path_save, "resolution_vs_centroid_genes.tsv"),
        sep="\t",
        index=False,
    )

    # Plot
    plt.figure(figsize=(6, 4))
    plt.plot(
        out["resolution"],
        out["n_centroid_genes"],
        marker="o",
    )
    plt.xlabel("Resolution")
    plt.ylabel("Number of centroid genes")
    plt.title("Model complexity vs resolution")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        Path(path_save, "resolution_vs_centroid_genes.png"),
        dpi=300,
    )
    plt.close()

def sankey_gene_centroid_flow(
    path_save,
    output_html="sankey_gene_centroid_flow.html",
    min_flow=1,
):
    """
    Build an interactive Sankey plot showing patient flow
    across gene-centroid clusters from high to low resolution.

    Nodes: GENE@RESOLUTION
    Links: number of shared patients between centroid clusters
           at consecutive resolutions.

    Output:
        sankey_gene_centroid_flow.html
    """

    import plotly.graph_objects as go
    import random
    import colorsys
    def _gene_color_map(genes):
        """
        Assign a stable color to each gene.
        """
        colors = {}
        for i, g in enumerate(sorted(genes)):
            hue = i / max(1, len(genes))
            rgb = colorsys.hsv_to_rgb(hue, 0.6, 0.85)
            colors[g] = f"rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})"
        return colors



    # --------------------------------------------------
    # Load centroid table
    # --------------------------------------------------
    df = pd.read_csv(
        Path(path_save, "resolution_centroids.tsv"),
        sep="\t",
    )

    # --------------------------------------------------
    # Parse clusters: resolution -> gene -> samples
    # --------------------------------------------------
    clusters = {}

    for col in df.columns:
        # Expected: gene_clusterX_res
        left, _, res = col.rpartition("_")
        gene_part = left.split("_cluster")[0]
        res = float(res)

        samples = set(
            s
            for s in df[col].dropna().astype(str)
            if s != ""
        )

        clusters.setdefault(res, {}).setdefault(
            gene_part, set()
        ).update(samples)

    # --------------------------------------------------
    # Sort resolutions (HIGH -> LOW)
    # --------------------------------------------------
    resolutions = sorted(clusters.keys(), reverse=True)

    # --------------------------------------------------
    # Assign colors to genes (constant across resolutions)
    # --------------------------------------------------
    all_genes = set()
    for res in clusters:
        all_genes.update(clusters[res].keys())

    gene_colors = _gene_color_map(all_genes)


    # --------------------------------------------------
    # Build nodes
    # --------------------------------------------------
    node_index = {}
    node_labels = []
    node_colors = []

    def _add_node(label):
        if label not in node_index:
            gene = label.split("@")[0]
            node_index[label] = len(node_labels)
            node_labels.append(label)
            node_colors.append(gene_colors.get(gene, "lightgrey"))


    # --------------------------------------------------
    # Build links
    # --------------------------------------------------
    sources = []
    targets = []
    values = []


    for r1, r2 in zip(resolutions[:-1], resolutions[1:]):
        for gene1, s1 in clusters[r1].items():
            for gene2, s2 in clusters[r2].items():
                overlap = len(s1 & s2)

                if overlap >= min_flow:
                    src = f"{gene1}@{r1}"
                    tgt = f"{gene2}@{r2}"

                    _add_node(src)
                    _add_node(tgt)

                    sources.append(node_index[src])
                    targets.append(node_index[tgt])
                    values.append(overlap)


    # --------------------------------------------------
    # Build Sankey plot
    # --------------------------------------------------
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=15,
                    line=dict(color="black", width=0.5),
                    label=node_labels,
                    color=node_colors,
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color="rgba(180,180,180,0.6)",
                ),
            )
        ]
    )

    fig.update_layout(
        title_text="Patient flow across gene-centric clusters (multi-resolution)",
        font_size=10,
    )

    # --------------------------------------------------
    # Save HTML
    # --------------------------------------------------
    out_path = Path(path_save, output_html)
    fig.write_html(out_path)

    print(f"[DONE] Sankey plot saved to {out_path}")
