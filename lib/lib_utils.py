import colorsys
import datetime
import os
import pickle
import random as random
import sys
from itertools import combinations
from multiprocessing import Pool

import igraph as ig
import numpy as np
import pandas as pd
import tap
from rpy2 import robjects


def load_df(config):
    """
    Load dataframes from CSV files based on the provided configuration.
    Args:
        config (dict): A dictionary containing file paths and parameters.
            The dictionary should have the following structure:
            {
                "Paths": {
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
        config["Paths"]["data_mutational"],
        sep=config["Paths"]["data_mutational_sep"],
        skiprows=config["Paths"]["data_mutational_skip"],
        low_memory=False,
    )
    data_clinical_sample = None
    if config["Paths"]["data_clinical_sample"] != "":
        data_clinical_sample = pd.read_csv(
            config["Paths"]["data_clinical_sample"],
            sep=config["Paths"]["data_clinical_sample_sep"],
            skiprows=config["Paths"]["data_clinical_sample_skip"],
            low_memory=False,
        )
    data_clinical_patient = None
    if config["Paths"]["data_clinical_patient"] != "":
        data_clinical_patient = pd.read_csv(
            config["Paths"]["data_clinical_patient"],
            sep=config["Paths"]["data_clinical_patient_sep"],
            skiprows=config["Paths"]["data_clinical_patient_skip"],
            low_memory=False,
        )
    return df_mut, data_clinical_sample, data_clinical_patient


def filter_vaf(config, df_mut):
    """
    Filters a DataFrame of mutations based on Variant Allele Frequency (VAF) score.
    This function filters the mutations in the DataFrame `df_mut` based on the VAF score
    specified in the `config` dictionary. If the `column_vaf` is not specified, it defaults
    to "t_AF" and calculates the VAF using the `libu.calculated_vaf` function.
    Args:
        config (dict): Configuration dictionary containing the following keys:
            - "Mutation": A dictionary with the keys:
                - "vaf_score": The VAF score threshold.
                - "vaf_column": The column name for the VAF score in the DataFrame.
            - "Paths": A dictionary with the keys:
                - "name_study": The name of the study.
                - "data_mutational_sep": The separator to use when saving the filtered data.
        df_mut (DataFrame): The DataFrame containing mutational data.
    Returns:
        DataFrame: The filtered DataFrame with mutations that meet the VAF score threshold.
    Saves:
        A CSV file of the filtered mutations to the path specified in the `config` dictionary.
    """
    column_vaf = config["Mutation"]["vaf_column"]

    if column_vaf != "":
        df_mut = df_mut[df_mut[column_vaf] >= config["Mutation"]["vaf_score"]]
    else:
        config["Mutation"]["vaf_column"] = "t_AF"
        df_mut["t_AF"] = df_mut.apply(calculated_vaf, axis=1)
        df_mut = df_mut[df_mut["t_AF"] >= config["Mutation"]["vaf_score"]]

    df_mut.to_csv(
        os.path.join(
            "study",
            config["Paths"]["name_study"],
            "input",
            "data_mutational_filtered.txt",
        ),
        index=False,
        sep=config["Paths"]["data_mutational_sep"],
    )
    return config, df_mut


def enrichment_with_r(path_save, map_cluster):
    """
    Perform enrichment analysis using an R script for given gene clusters.
    This function loads an R script and uses it to perform enrichment analysis
    on gene clusters. The results are saved in specified directories.
    Parameters:
    path_save (str): The path where the results will be saved.
    map_cluster (list): A list of cluster identifiers to be analyzed.
    Returns:
    None
    """

    # Load R script
    robjects.r.source("./lib/enrichment.r")
    # R functions
    r_func = robjects.globalenv["all_analisi"]
    g_path = os.path.join(path_save, "Gene")
    o_path = os.path.join(path_save, "Arricchimento_all_genes")
    os.makedirs(os.path.join(o_path, "GO"), exist_ok=True)
    os.makedirs(os.path.join(o_path, "KEGG"), exist_ok=True)
    os.makedirs(os.path.join(o_path, "WIKI"), exist_ok=True)
    os.makedirs(os.path.join(o_path, "REACTOME"), exist_ok=True)
    if sys.platform.startswith("win") or sys.platform.startswith("linux"):
        with Pool() as p:
            p.map(
                r_func,
                [
                    [os.path.join(g_path, f"genes_cluster_{c}.csv"), c, o_path]
                    for c in map_cluster
                ],
            )
    else:
        for c in map_cluster:
            r_func([os.path.join(g_path, f"genes_cluster_{c}.csv"), c, o_path])


def adding_category_mutation(data_mutational, list_columns):
    nuovi_nomi = {}
    # FILL NA
    for col in list_columns:
        data_mutational.fillna({f"{col}": "N/D"}, inplace=True)
    # GROUP
    gruppi_mutazioni = data_mutational.groupby(list_columns)
    # GENERATE ID
    for grp, data in gruppi_mutazioni:
        # GENE_CHROMOSOME_START_END
        nuovi_nomi[grp] = "Mut_" + "_".join(str(g) for g in grp)
    # ADD COLUMN
    data_mutational["TN_mutation_label"] = data_mutational.apply(
        lambda row: nuovi_nomi[tuple(row[c] for c in list_columns)], axis=1
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
    data_mutational, column_mutation, column_gene, column_sample, vaf, column_vaf
):
    map_variants = {}
    map_patients = {}
    for _index, _row in data_mutational.iterrows():
        # GET INFOS FROM DATA
        _paz = str(_row[column_sample])
        _gene = str(_row[column_gene])
        _category = str(_row[column_mutation])

        # se questa colonna è vuota tornerà una stringa vuota
        if "HGVSp_Short" in data_mutational.columns:
            _sost_amm = str(_row["HGVSp_Short"])
        else:
            _sost_amm = ""
        # se questa colonna è vuota tornerà una stringa vuota (dire all'utente che la colonna corrispondente alla VAF dovrà essere chiamata t_AF)
        if vaf:
            _vaf = float(_row[column_vaf])
        else:
            _vaf = ""

        # creazione dizionario delle varianti
        if _category not in map_variants.keys():
            map_variants[_category] = {}
            map_variants[_category]["patients"] = set()
            map_variants[_category]["sost_amm"] = _sost_amm
            map_variants[_category]["gene"] = _gene
            map_variants[_category]["vaf"] = _vaf

        map_variants[_category]["patients"].add(_paz)

        # creazione dizionario dei pazienti
        if _paz not in map_patients.keys():
            map_patients[_paz] = {}
            map_patients[_paz]["variants"] = set()
        map_patients[_paz]["variants"].add(_category)

    return map_patients, map_variants


# Creazione del Grafo
def graph_creation(map_patients, map_variants):
    edges = []
    for _k_variant, _v_variant in map_variants.items():
        for _k_patient, _v_patient in map_patients.items():
            if _k_patient in [_x for _x in _v_variant["patients"]]:
                edges.append((_k_variant, _k_patient))

    graph = ig.Graph()
    graph.add_vertices(
        list(map_variants.keys()),
        attributes={
            "vertex_type": "VARIANT",
            "color_vertex": "blue",
            "shape_vertex": "circle",
            "gene": [f'{value["gene"]}' for key, value in map_variants.items()],
            "sost_amm": [f'{value["sost_amm"]}' for key, value in map_variants.items()],
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
            if vertex["gene"] not in gene_total_count.keys():
                gene_total_count[vertex["gene"]] = 0
            gene_total_count[vertex["gene"]] += 1
    gene_total_count = dict(
        sorted(gene_total_count.items(), key=lambda kv: kv[1], reverse=True)
    )
    return gene_total_count


def process_data(args):
    GRAPH = args[0]
    _seed = args[1]
    random.seed(_seed)
    _dendro_2 = GRAPH.community_leiden(objective_function="modularity")
    modularity = _dendro_2.modularity
    return modularity


# SELEZIONE DEL SEED CHE Dà VALORE DI MODULARITà PIù ALTA A SEGUITO DEL LEIDEN ALGORITHM
def selected_seed(GRAPH):
    best_seed = 0
    if sys.platform.startswith("win") or sys.platform.startswith("linux"):
        data = []
        for s in range(1000):
            data.append((GRAPH, s))
        with Pool() as p:
            mod_results = p.map(process_data, data)
        best_seed = mod_results.index(max(mod_results))
    else:
        mod_results = []
        data = []
        for s in range(1000):
            data.append((GRAPH, s))
        for s in data:
            mod_results.append(process_data(s))
        best_seed = mod_results.index(max(mod_results))

    best_seed = mod_results.index(max(mod_results))
    return best_seed


# SELEZIONE DEL SEED E LANCIO DELL'ALGORITMO DI CLUSTERIZZAZIONE
def leiden_clustering(graph, best_seed):
    random.seed(best_seed)
    dendro = graph.community_leiden(objective_function="modularity")
    print("numero di clusters:", len(list(dendro)), "Modularità:", dendro.modularity)
    return dendro


# adding color for cluster (for cytoscape)
def adding_graph_color(graph, dendro):
    num_clusters = len(dendro)
    # Genera colori in base al numero di cluster
    colors = []
    for i in range(num_clusters):
        hue = i / num_clusters  # variare la tonalità in base al numero di cluster
        rgb = colorsys.hsv_to_rgb(hue, 1, 1)  # converti da spazio colore HSV a RGB
        colors.append(
            "#%02x%02x%02x" % tuple(int(c * 255) for c in rgb)
        )  # formato esadecimale RGB
    # Assegna i colori alle singole comunità nel grafo
    graph.vs["color"] = [colors[cluster] for cluster in dendro.membership]

    return graph


# FIXME
def plot_graph(graph, path_save, gene):
    ig.plot(
        graph,
        f"{path_save}/plot_{gene}.pdf",
        **{
            "vertex_color": graph.vs["color"],
            "bbox": (2000, 1000),
            "edge_curved": 0.1,
            "vertex_label": graph.vs["name"],
            "vertex_label_size": 1.5,
            "vertex_size": 20,
            "edge_color": "grey",
            "vertex_shape": [
                "triangle" if v["vertex_type"] == "PATIENT" else "circle"
                for v in graph.vs
            ],
            # "edge_color": [get_color(x, _max_conn, _min_conn) for x in _new_g.es["weight"]],
            "edge_widht": 0.3,
            "vertex_frame_width": 0.05,
            "layout": "fr",
        },
    )


# function to plot single graph
def plot_graph_single_graph(
    g,
    path_save,
    label="name",
    color="color",
    shape="vertex_shape",
    layout="kk",
    title="plot_graph",
    on_file=True,
):
    visual_style = {}
    visual_style["layout"] = g.layout(layout)
    visual_style["vertex_label"] = g.vs[label]
    visual_style["color"] = g.vs[color]
    visual_style["vertex_shape"] = g.vs["vertex_shape"]
    visual_style["vertex_label_size"] = 8
    visual_style["vertex_size"] = 35
    visual_style["edge_width"] = 0.01
    visual_style["edge_color"] = "#e3e0de"
    visual_style["bbox"] = (len(g.vs) * 25, len(g.vs) * 18)
    visual_style["margin"] = 45
    visual_style["color"] = g.vs[color]

    if on_file:
        os.makedirs(f"{path_save}/Images_Graph_Plot/", exist_ok=True)
        ig.plot(
            g,
            f"{path_save}/Images_Graph_Plot/{title}_{datetime.datetime.now().strftime('%Y_%m_%d')}.pdf",
            **visual_style,
        )
    else:
        ig.plot(g, **visual_style)


# function to add single vertex
def add_unique_vertex(g, value, kwds={}, debug=False):
    try:
        if g.vs.find(name=value):
            if debug:
                print(f"Vertex {value} already present")
            return
    except:
        # id not present
        pass

    g.add_vertex(value, **kwds)


# function to add single edge
def add_unique_edge(g, id1, id2, kwds={}, direct=False, debug=False):
    for e in g.es:
        if (
            g.vs[e.source]["name"] == id1
            and g.vs[e.target]["name"] == id2
            or not direct
            and g.vs[e.source]["name"] == id2
            and g.vs[e.target]["name"] == id1
        ):
            if debug:
                print("edge is already present")
            return
    try:
        g.add_edge(id1, id2)
    except:
        if debug:
            print("vertex id not valid")


# function to plot single cluster as a graph
def plot_cluster_as_graph(g, cluster_index, path_save):
    g_cluster = ig.Graph()
    for v in g.vs:
        if v["cluster"] == cluster_index:
            add_unique_vertex(
                g_cluster,
                v["name"],
                {"color": v["color_vertex"], "vertex_shape": v["shape_vertex"]},
            )
    for e in g.es:
        patient_cluster = g.vs[e.source]
        variant_cluster = g.vs[e.target]
        if (
            patient_cluster["cluster"] == cluster_index
            and variant_cluster["cluster"] == cluster_index
        ):
            add_unique_edge(g_cluster, patient_cluster["name"], variant_cluster["name"])

    plot_graph_single_graph(
        g_cluster, path_save, layout="kk", title=f"graph_cluster_{cluster_index}"
    )


# function to save graph as graphml file for cytoscape
def save_graph_to_file(graph, path_save):
    graph.write_graphml(f"{path_save}/grafo_cytoscape.graphml")
    with open(f"{path_save}/graph.pickle", "wb") as f:
        pickle.dump(graph, f, protocol=pickle.HIGHEST_PROTOCOL)


# function to create a map for cluster
def map_cluster_creation(graph, dendro):
    map_cluster = {}
    for cluster_index in range(len(dendro)):
        map_cluster[cluster_index] = []
        for element in dendro[cluster_index]:
            _vertex = graph.vs()[element]
            map_cluster[cluster_index].append((_vertex["name"], _vertex["vertex_type"]))
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
def centroids_cluster(dendro, path_save):
    with open(f"{path_save}/Centroidi_Mutazioni.csv", "w") as f:
        for _i in range(len(dendro)):
            sub_graph = dendro.subgraph(_i)
            max_value = 0
            _list_centroids = []
            for _v in sub_graph.vs():
                # print(_v["name"])
                if _v["vertex_type"] != "VARIANT":
                    continue
                temp_val = sub_graph.neighborhood_size(_v)
                if temp_val > max_value:
                    max_value = temp_val
                    _list_centroids = [(_v["name"], max_value - 1)]
                elif temp_val == max_value:
                    _list_centroids.append((_v["name"], max_value - 1))
            # print(path_save)
            f.write(
                f"Cluster {_i} centroids found {len(_list_centroids)}: {_list_centroids}\n"
            )


# funzione per scrivere e trovare le coppie centroide-elemento per ciascun cluster
def couple_centroid_element(dendro, map_cluster, path_save):
    all_pairs = []
    for _i in range(len(dendro)):
        sub_graph = dendro.subgraph(_i)
        max_value = 0
        _list_centroids = []
        # individuazione centroide per ogni cluster
        for _v in sub_graph.vs():
            temp_val = sub_graph.neighborhood_size(_v)
            if temp_val > max_value:
                max_value = temp_val
                _list_centroids = [(_v["name"], max_value - 1)]
            elif temp_val == max_value:
                _list_centroids.append((_v["name"], max_value - 1))

        _list_pairs = []
        # creazione coppie centroide-elemento per ogni cluster
        for element in map_cluster[_i]:
            if element[0] not in [x[0] for x in _list_centroids]:
                for c in _list_centroids:
                    _list_pairs.append((c[0], element[0]))
        all_pairs += _list_pairs
    with open(f"{path_save}/coppie_centroide.csv", "w") as f:
        f.write("Centroide\tElemento\n")
        for elements in all_pairs:
            f.write(elements[0] + "\t" + elements[1] + "\n")
    return all_pairs


# funzione per scrivere il numero di connessioni che ciascuna variante ha nel cluster
def degree_variant_cluster(map_cluster, graph, path_save):
    os.makedirs(f"{path_save}/Variants_Degree", exist_ok=True)
    for cluster_index in map_cluster.keys():
        list_vertices_filtered = graph.vs.select(
            lambda x: x["cluster"] == cluster_index
        )
        g_cluster = graph.induced_subgraph(list_vertices_filtered)
        degrees = g_cluster.degree()
        with open(
            f"{path_save}/Variants_Degree/variants_degree_cluster{cluster_index}.csv",
            "w",
        ) as f:
            f.write("Variants\tDegree\n")
            for i, degree in enumerate(degrees):
                if g_cluster.vs[i]["vertex_type"] != "PATIENT":
                    f.write(f"{g_cluster.vs[i]['name']}\t{degree}\n")


# AGGIUNTA DELLE INFORMAZIONI CLINICHE DI INTERESSE ALLA MAPPA DEI PAZIENTI + AGGIUNTA DEL CLUSTER DI APPARTENENZA ALLA MAPPA DEI PAZIENTI E DELLE MITAZIONI
def enriched_sample_data(data_sample, MAP_PATIENTS, sample_name):
    # generalizzazione --> CLINICAL SAMPLE
    list_clinical_parameters = [col for col in data_sample.columns]
    for _i, _row in data_sample.iterrows():
        _sample = _row[sample_name]
        for parameter in list_clinical_parameters:
            variable = _row[parameter]
            if _sample in MAP_PATIENTS.keys():
                MAP_PATIENTS[_sample][parameter] = variable
    return MAP_PATIENTS


def enriched_patient_data(data_patient, MAP_PATIENTS, patient_name):
    # generalizzazione --> CLINICAL PATIENT
    list_clinical_patient = [col for col in data_patient.columns]
    for _i, _row in data_patient.iterrows():
        _paz = _row[patient_name]
        for parameters in list_clinical_patient:
            variable = _row[parameters]
            for _sample, value in MAP_PATIENTS.items():
                if value[patient_name] == _paz:
                    MAP_PATIENTS[_sample][parameters] = variable
    return MAP_PATIENTS


# function to add clinical information to the vertex
def adding_clinical_info_graph(graph, map_patients):
    for vertex in graph.vs:
        nome_paziente = vertex["name"]
        if nome_paziente in map_patients.keys():
            # Aggiungi le info cliniche come attributi del nodo
            for key, value in map_patients[nome_paziente].items():
                vertex[key] = value
    return graph


# CREAZIONE DEL FILE "CLUSTER_CLINICAL_DATA" IN CUI INSERIAMO IL CLUSTER DI APPARTENENZA DI OGI SAMPLES + INFORMAZIONI CLINICHE
def creation_cluster_clinical_data(MAP_PATIENTS, path_saved):
    header_clinical_data = sorted(
        list(
            set(
                [
                    k2
                    for k, v in MAP_PATIENTS.items()
                    for k2 in v.keys()
                    if k2 != "variants"
                ]
            )
        )
    )
    with open(f"{path_saved}/cluster_clinical_data.csv", "w") as f:
        f.write("\t".join(header_clinical_data) + "\n")
        for k, v in MAP_PATIENTS.items():
            temp = []
            for data in header_clinical_data:
                temp.append(str(v[data]))
            f.write("\t".join(temp) + "\n")


# creazione di un file con il numero di pazienti e varianti per ogni cluster  + CREAZIONE DI UN ARRAY CON I CLUSTER CON > 1 PAZIENTE + UN ARRAY CON I CLUSTER CON ==1 PAZIENTE
def cluster_division(MAP_CLUSTER):
    cluster_one_patient = []
    cluster_more_patient = []
    with open("numerosity_cluster.csv", "w") as f:
        f.write("Cluster\tPatient\tVariant\n")
        for cluster, value in MAP_CLUSTER.items():
            count_variant = 0
            count_patient = 0
            for v in value:
                if v[1] == "VARIANT":
                    count_variant += 1
                else:
                    count_patient += 1
            f.write(
                str(cluster)
                + "\t"
                + str(count_patient)
                + "\t"
                + str(count_variant)
                + "\n"
            )
            if count_patient == 1:
                cluster_one_patient.append(cluster)
            else:
                cluster_more_patient.append(cluster)

    return cluster_more_patient, cluster_one_patient


# defizione del numero di connessioni per ogni variante all'interno del cluster
def variant_conncection_patient(_dendro_2):
    variant_patient_connection_count = {}
    for _i in range(len(_dendro_2)):
        variant_patient_connection_count[_i] = []
        sub_graph = _dendro_2.subgraph(_i)
        for _v in sub_graph.vs():
            if _v["vertex_type"] != "VARIANT":
                continue
            temp_val = sub_graph.neighborhood_size(_v)
            variant_patient_connection_count[_i].append((_v["name"], (temp_val - 1)))
    variant_patient_connection_count = {
        key: sorted(value, key=lambda x: x[1], reverse=True)
        for key, value in variant_patient_connection_count.items()
    }
    return variant_patient_connection_count


# defizione del numero di connessioni per ogni paziente all'interno del cluster
def patient_connection_variant(_dendro_2):
    patient_variant_connection_count = {}
    for _i in range(len(_dendro_2)):
        patient_variant_connection_count[_i] = []
        sub_graph = _dendro_2.subgraph(_i)
        for _v in sub_graph.vs():
            if _v["vertex_type"] != "PATIENT":
                continue
            temp_val = sub_graph.neighborhood_size(_v)
            patient_variant_connection_count[_i].append((_v["name"], (temp_val - 1)))
    patient_variant_connection_count = {
        key: sorted(value, key=lambda x: x[1], reverse=True)
        for key, value in patient_variant_connection_count.items()
    }
    return patient_variant_connection_count


# CREAZIONE DI UN FILE "CONNECTION_VARIANT" IN CUI SONO INDICATE IL NUMERO DI VARIANTI COMUNI TRA I VARI PAZIENTI DI UN CLUSTER
def file_connection_variant(MAP_CLUSTER, MAP_PATIENTS, path_saved):
    with open(f"{path_saved}/connection_variant.csv", "w") as f:
        f.write("Cluster\tPatient_1\tPatient_2\tVariant\tNumber_variant\n")
        for cluster, values in MAP_CLUSTER.items():
            # f.write(str(cluster)+"\t")
            patients = []
            for value in values:
                if value[1] == "PATIENT":
                    patients.append(value[0])
            for patient in list(combinations(patients, 2)):
                paz_1 = patient[0]
                paz_2 = patient[1]
                variant_p1 = MAP_PATIENTS[paz_1]["variants"]
                variant_p2 = MAP_PATIENTS[paz_2]["variants"]
                variant_common = list(variant_p1.intersection(variant_p2))
                f.write(
                    str(cluster)
                    + "\t"
                    + str(paz_1)
                    + "\t"
                    + str(paz_2)
                    + "\t"
                    + ",".join(variant_common)
                    + "\t"
                    + str(len(variant_common))
                    + "\n"
                )


# CREAZIONE DI UN FILE "CONCCECTION_PATIENT" IN CUI SONO INDICATI IL NUMERO DI PAZIENTI COMUNI TRA LE VARIE VARIANTI DI UN CLUSTER
def file_connection_patient(MAP_CLUSTER, MAP_VARIANTS, path_saved):
    with open(f"{path_saved}/connection_patient.csv", "w") as f:
        f.write("Cluster\tVariant_1\tVariant_2\tPatient\tNumber_patient\n")
        for cluster, values in MAP_CLUSTER.items():
            # f.write(str(cluster)+"\t")
            variants = []
            for value in values:
                if value[1] == "VARIANT":
                    variants.append(value[0])
            for variant in list(combinations(variants, 2)):
                var_1 = variant[0]
                var_2 = variant[1]
                patient_v1 = MAP_VARIANTS[var_1]["patients"]
                patient_v2 = MAP_VARIANTS[var_2]["patients"]
                patient_common = list(patient_v1.intersection(patient_v2))
                f.write(
                    str(cluster)
                    + "\t"
                    + str(var_1)
                    + "\t"
                    + str(var_2)
                    + "\t"
                    + ",".join(patient_common)
                    + "\t"
                    + str(len(patient_common))
                    + "\n"
                )


# funzione che scrive un file che riassume le info presenti in ogni cluster
def numerosity_info(g, map_cluster, path_save):
    with open(f"{path_save}/numerosity_cluster.csv", "w") as f:
        f.write("Cluster\tPatient\tVariant\tGene\n")
        for cluster, value in map_cluster.items():
            count_variant = len(
                [
                    v
                    for v in g.vs
                    if v["cluster"] == cluster and v["vertex_type"] == "VARIANT"
                ]
            )
            count_patient = len(
                [
                    v["gene"]
                    for v in g.vs
                    if v["cluster"] == cluster and v["vertex_type"] == "PATIENT"
                ]
            )
            gene_count = len(
                set(
                    v["gene"]
                    for v in g.vs
                    if v["cluster"] == cluster and v["vertex_type"] == "VARIANT"
                )
            )
            f.write(
                str(cluster)
                + "\t"
                + str(count_patient)
                + "\t"
                + str(count_variant)
                + "\t"
                + str(gene_count)
                + "\n"
            )


# riassunto delle informazioni (mutazioni e pazienti per ciascun cluster)
def summary_info(g, map_cluster, patient_column, path_save):
    with open(f"{path_save}/summury_file.csv", "w") as f:
        f.write("Cluster\tNode_Name\tType\tGene\tPatient_id\n")
        for cluster, value in map_cluster.items():
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
                        + "\n"
                    )
                else:
                    f.write(
                        str(cluster)
                        + "\t"
                        + str(v["name"])
                        + "\tPATIENT\t"
                        + "None"
                        + "\t"
                        + str(v[patient_column])
                        + "\n"
                    )


# aggiunta size ai vertici mutazione
def add_size_node(GRAPH, variant_patient_connection_count):
    for cluster, elements in variant_patient_connection_count.items():
        for e in elements:
            variant = e[0]
            presence = e[1]
            _vx = GRAPH.vs.find(name=variant)
            _vx["size"] = presence
    return GRAPH


# COUNT DEI GENI PRESENTI IN OGNI CLUSTER + VALORE PERCENTUALE DI APPARTENENZA DI OGNI GENE A OGNI CLUSTER, CON SALVATAGGIO DELLE INFO IN UN FILE "DISTRIBUTION_GENE_CLUSTER"
def count_gene_abs_percent(g, map_cluster, gene_total_count, path_save):
    MAP_CLUSTER_GENE = {}
    for cluster, values in map_cluster.items():
        MAP_CLUSTER_GENE[cluster] = {}
        _list_genes = [
            v["gene"]
            for v in g.vs
            if v["cluster"] == cluster and v["vertex_type"] == "VARIANT"
        ]
        for _g in set(_list_genes):
            count = _list_genes.count(_g)
            MAP_CLUSTER_GENE[cluster][_g] = count

    MAP_CLUSTER_GENE_PERCEN = {}
    for cluster, genes in MAP_CLUSTER_GENE.items():
        MAP_CLUSTER_GENE_PERCEN[cluster] = {}
        for gene, count in genes.items():
            if gene not in MAP_CLUSTER_GENE_PERCEN[cluster].keys():
                MAP_CLUSTER_GENE_PERCEN[cluster][gene] = (
                    MAP_CLUSTER_GENE[cluster][gene] * 100
                ) / gene_total_count[gene]

    MAP_CLUSTER_GENE_PERCEN = {
        k: dict(sorted(v.items(), key=lambda x: x[1], reverse=True))
        for k, v in MAP_CLUSTER_GENE_PERCEN.items()
    }
    with open(f"{path_save}/distribution_gene_cluster.csv", "w") as f:
        f.write("Cluster\tGene\tPercentage\n")
        for cluster, genes in MAP_CLUSTER_GENE_PERCEN.items():
            for gene, count in genes.items():
                f.write(str(cluster) + "\t" + str(gene) + "\t" + str(count) + "\n")

    return MAP_CLUSTER_GENE, MAP_CLUSTER_GENE_PERCEN


# creazione di un file per ogni cluster, contenente i geni presenti
def genes_single_cluster(g, map_cluster, path_save):
    os.makedirs(f"{path_save}/Gene", exist_ok=True)
    for cluster, infos in map_cluster.items():
        with open(f"{path_save}/Gene/genes_cluster_{cluster}.csv", "w") as f:
            set_gene = set(
                v["gene"]
                for v in g.vs
                if v["cluster"] == cluster and v["vertex_type"] == "VARIANT"
            )
            for genes in set_gene:
                f.write(genes + "\n")


# creazione di un file che per ogni cluster, tiene il conto del numero di mutazioni presenti sul gene
def genes_count_mutation_single_cluster(map_cluster_gene_abs, path_save):
    os.makedirs(f"{path_save}/Gene_Count", exist_ok=True)
    for cluster, infos in map_cluster_gene_abs.items():
        with open(f"{path_save}/Gene_Count/genes_cluster_{cluster}.csv", "w") as f:
            f.write("GENE\tCOUNT\n")
            for k, v in infos.items():
                f.write(f"{k}\t{v}\n")


# AGGIUNTA DELL'ATTRIBUTO CLUSTER AI NODI DEL GRAFO
def cluster_noded_attributes(GRAPH, MAP_PATIENTS, MAP_VARIANTS):
    for v in GRAPH.vs():
        if v["vertex_type"] == "PATIENT":
            v["cluster"] = MAP_PATIENTS[v["name"]]["cluster"]
        else:
            v["cluster"] = MAP_VARIANTS[v["name"]]["cluster"]
    return GRAPH


# CREAZIONE, PER OGNI CLUSTER, DI UN GRAFO I CUI NODI SONO LE VARIANTI E IN CUI:
# - LA GRANDEZZA DEI NODI è = AL NUMERO DI PAZIENTI DEL CLUSTER CHE HANNO QUELLA MUTAZIONE
# - 2 MUTAZIONI SONO CONNESSE SE CO-MUTATE (QUINDI PRESENTI IN ALMENO 2 PAZIENTI DEL CLUSTER)
# - GRANDEZZA DEI NODI è PARI AL NUMERO DI PAZIENTI CHE HANNO QUELLE DUE MUTAZIONI NEL CLUSTER
def get_color_comutated_cluster(value, _max, _min, _scaling=1):
    _range = (_max + 1 - _min) * _scaling
    _range_sector = _range / 4

    _index_color = (value + 1) // _range_sector

    if _index_color == 0:
        # RED
        return (1.0, 0.0, 0.0, 0.5)
    elif _index_color == 1:
        # ORANGE
        return (1.0, 0.5, 0.0, 0.5)
    elif _index_color == 2:
        # YELLOW
        return (1.0, 1.0, 0.0, 0.5)
    else:
        # GREEN
        return (0.0, 1.0, 0.0, 0.5)


def plot_distance_comutated_cluster_variants(
    _dendro_2, cluster_index, connection_df, path_saved
):
    os.makedirs(f"{path_saved}/plot_mutation_connection", exist_ok=True)

    sub_graph = _dendro_2.subgraph(cluster_index)
    list_variants = list(sub_graph.vs.select(vertex_type="VARIANT"))
    _new_g = ig.Graph()
    for var in list_variants:
        _new_g.add_vertex(
            name=var["name"],
            **{"size": min(100, var["size"] * 5), "sost_amm": var["sost_amm"]},
        )
    _temp_edge = set()

    _df_cluster = connection_df[connection_df["Cluster"] == cluster_index]
    _max_conn = max(_df_cluster["Number_patient"].values)
    _min_conn = min(_df_cluster["Number_patient"].values)

    for var in list_variants:
        _var_name = var["name"]
        # _var_name = var["sost_amm"]
        _filtered_edge = _df_cluster[
            (_df_cluster["Variant_1"] == _var_name)
            | (_df_cluster["Variant_2"] == _var_name)
        ]

        for index, row in _filtered_edge.iterrows():
            _var1 = row["Variant_1"]
            _var2 = row["Variant_2"]
            _w = row["Number_patient"]
            if _w > 0:
                if _var1 > _var2:
                    _temp_edge.add((_var1, _var2, _w))
                else:
                    _temp_edge.add((_var2, _var1, _w))

    for e in _temp_edge:
        _new_g.add_edge(e[0], e[1], weight=e[2])

    weight_arch = []

    if len(_new_g.es) == 0:
        weight_arch = []
    else:
        weight_arch = _new_g.es["weight"]

    ig.plot(
        _new_g,
        # target=ax,
        f"{path_saved}/plot_mutation_connection/cluster_{cluster_index}.pdf",
        **{
            # "edge_width": _new_g.es["weight"],
            "edge_width": weight_arch,
            "vertex_color": "cyan",
            "bbox": (1290, 820),
            "edge_curved": 0.1,
            # "vertex_label": [x.replace("_", "\n") for x in _new_g.vs["name"]],
            "vertex_label": [x.replace("_", "\n") for x in _new_g.vs["sost_amm"]],
            "vertex_label_size": 2.5,
            # "edge_color": [get_color(x, _max_conn, _min_conn) for x in _new_g.es["weight"]],
            "edge_color": [
                get_color_comutated_cluster(x, _max_conn, _min_conn)
                for x in weight_arch
            ],
            "edge_label": weight_arch,
            "edge_label_size": [1.5 + (n * 2) for n in weight_arch],
            "vertex_frame_width": 0.05,
            "background": (0.3, 0.4, 0.5, 1),
        },
    )


def plot_distance_comutated_cluster_patients(
    _dendro_2, cluster_index, connection_df, path_saved
):
    os.makedirs(f"{path_saved}/plot_patient_connection", exist_ok=True)

    sub_graph = _dendro_2.subgraph(cluster_index)
    list_patients = list(sub_graph.vs.select(vertex_type="PATIENT"))
    _new_g = ig.Graph()
    for paz in list_patients:
        _new_g.add_vertex(name=paz["name"])
    _temp_edge = set()

    _df_cluster = connection_df[connection_df["Cluster"] == cluster_index]
    _max_conn = max(_df_cluster["Number_variant"].values)
    _min_conn = min(_df_cluster["Number_variant"].values)

    for paz in list_patients:
        _paz_name = paz["name"]
        # _var_name = var["sost_amm"]
        _filtered_edge = _df_cluster[
            (_df_cluster["Patient_1"] == _paz_name)
            | (_df_cluster["Patient_2"] == _paz_name)
        ]

        for index, row in _filtered_edge.iterrows():
            _var1 = row["Patient_1"]
            _var2 = row["Patient_2"]
            _w = row["Number_variant"]
            if _w > 0:
                if _var1 > _var2:
                    _temp_edge.add((_var1, _var2, _w))
                else:
                    _temp_edge.add((_var2, _var1, _w))

    for e in _temp_edge:
        _new_g.add_edge(e[0], e[1], weight=e[2])

    weight_arch = []

    if len(_new_g.es) == 0:
        weight_arch = []
    else:
        weight_arch = _new_g.es["weight"]

    ig.plot(
        _new_g,
        # target=ax,
        f"{path_saved}/plot_patient_connection/cluster_{cluster_index}.pdf",
        **{
            # "edge_width": _new_g.es["weight"],
            "edge_width": weight_arch,
            "vertex_color": "cyan",
            "bbox": (1280, 720),
            "edge_curved": 0.1,
            # "vertex_label": [x.replace("_", "\n") for x in _new_g.vs["name"]],
            "vertex_label": _new_g.vs["name"],
            "vertex_label_size": 4,
            # "edge_color": [get_color(x, _max_conn, _min_conn) for x in _new_g.es["weight"]],
            "edge_color": [
                get_color_comutated_cluster(x, _max_conn, _min_conn)
                for x in weight_arch
            ],
            "edge_label": weight_arch,
            "edge_label_size": [5 + (n * 2) for n in weight_arch],
            "vertex_frame_width": 0.05,
            "background": (0.3, 0.4, 0.5, 1),
        },
    )


# prendi paziente 0 del cluster A
# calcolo indice di somiglianza con tutti i glia altri pazienti contenuti nel cluster B (tranne se stesso)
# i valori vengono aggiunti ad un array temporaneo
# passo al paziente successivo, una volta terminato il cliclo avrò un array con i valori di somiglianza
# return


def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(set(list1)) + len(set(list2))) - intersection
    return float(intersection) / union


def somiglianze_cluster(cluster_a, cluster_b, MAP_CLUSTER, MAP_PATIENTS):
    list_patients_cluster_a = [x for x in MAP_CLUSTER[cluster_a] if x[1] == "PATIENT"]
    list_patients_cluster_b = [x for x in MAP_CLUSTER[cluster_b] if x[1] == "PATIENT"]
    values = []

    # IN
    if cluster_a == cluster_b:
        for paz in list(combinations(list_patients_cluster_a, 2)):
            variant_paz_a = MAP_PATIENTS[paz[0][0]]["variants"]
            variant_paz_b = MAP_PATIENTS[paz[1][0]]["variants"]
            values.append(jaccard_similarity(variant_paz_a, variant_paz_b))
    # OUT
    else:
        for paz_a in list_patients_cluster_a:
            for paz_b in list_patients_cluster_b:
                variant_paz_a = MAP_PATIENTS[paz_a[0]]["variants"]
                variant_paz_b = MAP_PATIENTS[paz_b[0]]["variants"]
                values.append(jaccard_similarity(variant_paz_a, variant_paz_b))

    return values


def box_plot_similitudine_one_cluster(cluster_more_patient, path_saved):
    cluster_base = 0
    values = somiglianze_cluster(cluster_base, cluster_base)
    values_out = []
    for i in cluster_more_patient:
        if i == cluster_base:
            continue
        values_out += somiglianze_cluster(cluster_base, i)
    df_temp = pd.DataFrame(
        {
            "cluster_base": [cluster_base] * len(values)
            + [cluster_base] * len(values_out),
            "category": [f"C_{cluster_base}_IN"] * len(values)
            + [f"C_{cluster_base}_OUT"] * len(values_out),
            "values": values + values_out,
        }
    )
    print(values)
    tap.plot_stats(
        df_temp,
        "category",
        "values",
        filename=f"{path_saved}/adhesion_cluster_{cluster_base}.png",
    )


def box_plot_similitudine_all_clusters(
    cluster_more_patient,
):
    for index_cluster in cluster_more_patient:
        # IN
        df_values = pd.DataFrame(columns=["Cluster", "Category", "Value"])
        values_in = somiglianze_cluster(index_cluster, index_cluster)
        for value in values_in:
            df_values.loc[len(df_values.index)] = [
                index_cluster,
                f"C_{index_cluster}_IN",
                value,
            ]
        # OUT
        for index_cluster_out in cluster_more_patient:
            if index_cluster == index_cluster_out:
                continue
            values_out = somiglianze_cluster(index_cluster, index_cluster_out)
            for value in values_out:
                df_values.loc[len(df_values.index)] = [
                    index_cluster,
                    f"C_{index_cluster}_OUT",
                    value,
                ]

        tap.plot_stats(df_values, "Category", "Value")
