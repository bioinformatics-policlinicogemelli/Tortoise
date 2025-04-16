#!/usr/bin/env python3
"""Dash web application for visualizing and analyzing data.

The application includes multiple pages for different types of analysis,
including study creation, study description, pathway analysis,
clinical data visualization, survival analysis, and cluster comparison.

Modules:
    - base64
    - io
    - json
    - os
    - dash_bootstrap_components as dbc
    - dash_cytoscape as cyto
    - matplotlib
    - matplotlib.pyplot as plt
    - numpy as np
    - pandas as pd
    - plotly.express as px
    - plotly.graph_objects as go
    - tap
    - dash
    - lifelines
    - lifelines.statistics
    - plotly.subplots
    - lib.venn as venn
Functions:
    - filter_graph(cluster): Filters the graph based on the selected cluster.
    - redirect_pages(pathname): Redirects to the page based on the URL.
    - update_dropdown_liststudy(): Updates the dropdown list of studies.
    - select_study(value): Selects a study and loads its data.
    - update_mutational_file(filename): Updates mutational file component.
    - update_clinical_patient_file(filename): Updates patient file component.
    - update_clinical_sample_file(filename): Updates sample file component.
    - update_val_col_patient_name(data, sep, skip): Updates dropdown options.
    - update_val_col_sample_name(data, sep, skip): Updates dropdown options.
    - update_list_columns_mutation(data, sep, skip): Updates dropdown options.
    - create_study(...): Creates a new study based on the provided data.
    - dropdpwn_cluster(): Updates the cluster dropdown options.
    - display_node_data(data_dict): Displays data for the selected node.
    - update_graph(layout): Updates the graph layout.
    - update_cluster(cluster): Updates the cluster elements and figures.
    - update_go(cluster, pvalue, adj_pvalue, p_type): Updates the GO plot.
    - update_kegg(cluster, pvalue, adj_pvalue): Updates KEGG plot.
    - update_reactome(cluster, pvalue, adj_pvalue): Updates REACTOME plot.
    - update_wiki(cluster, pvalue, adj_pvalue): Updates WIKI plot.
    - dropdown_box_fig_1(): Updates dropdown options.
    - dropdown_box_fig_2(): Updates dropdown options.
    - func_single_plot(cluster, col_name): Generates a single cluster plot.
    - update_box_1(cluster, col_name): Updates the first box plot.
    - update_box_2(cluster, col_name): Updates the second box plot.
    - update_table_clinical_data(cluster): Updates the clinical data table.
    - update_overall_survival(cluster): Updates survival plot and statistics.
    - update_survival_comparison(list_clusters): Updates the survival plot.
    - update_venn(list_clusters): Updates the Venn diagram for gene.
    - update_genes_common(list_clusters): Updates the table of common genes.
    - func_multi_plot(list_clusters, col_name): Generates multi-cluster plot.
    - update_multi_fig(list_clusters, col_name): Updates multi-cluster plot.
Usage:
    Run this script to start the Dash web application.
    The application will be available at http://127.0.0.1:8593.
"""

import base64
import io
import json
import os
from pathlib import Path

import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_uploader as du
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import tap
from dash import (
    Dash,
    Input,
    Output,
    State,
    callback,
    dash_table,
    dcc,
    html,
    no_update,
)
from lifelines import KaplanMeierFitter
from lifelines.statistics import pairwise_logrank_test
from plotly.subplots import make_subplots
from rpy2.robjects import conversion, default_converter

from lib import venn

# Config lib
mpl.use("agg")
pd.options.mode.copy_on_write = True
cyto.load_extra_layouts()
# Global Vars
APP_NAME = "TORTOISE"
APP_LOGO = "/assets/tortoise.png"
DF_CLINICAL_DATA = None
NUMERIC_COLUMNS_CLINICAL = []
ALL_COLUMNS_CLINICAL = []
GRAPH = None
CLUSTERS_INDEX = []
CLUSTER_SELECTED = None
CLUSTER_SELECTED_MULTI = []
BOX_FIG_SELECTED_1 = None
BOX_FIG_SELECTED_2 = None
PATH_CONFIG = None
READY_FOR_CREATION = False
CONTEXT_DATA = {
    "config": {},
    "name_study": None,
    "name_study_input": None,
    "list_studies": [],
    "out_root_path": None,
    "stats": {
        "num_patient": 0,
        "num_gene": 0,
        "num_variant": 0,
        "num_cluster": 0,
        "modularity": 0,
    },
}


def filter_graph(cluster):
    if GRAPH is None:
        return []
    list_vertices = GRAPH.vs.select(lambda x: x["cluster"] == cluster)
    graph_filtered = GRAPH.induced_subgraph(list_vertices)
    g_ele = []
    # convert and add vertex
    for e in graph_filtered.vs():
        _map = e.attributes()
        _map["id"] = e.index
        _map["variants"] = None
        g_ele.append(
            {
                "data": _map,
                "classes": e["vertex_type"],
                "id": e.index,
                "grabbable": False,
            },
        )
    # convert and add edges
    g_ele.extend(
        [
            {"data": {"source": e.source, "target": e.target}}
            for e in graph_filtered.es()
        ],
    )
    return g_ele


# APP + SIDEBAR
APP = Dash(
    title=APP_NAME,
    external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME],
)
du.configure_upload(APP, "temp/", use_upload_id=False)
# ICON -> https://fontawesome.com/search
SIDEBAR = html.Div(
    [
        html.Div(
            [
                html.Img(src=APP_LOGO, className="navbar_logo"),
                html.Span(APP_NAME, className="navbar_title"),
            ],
            className="sidebar-header",
        ),
        html.Hr(),
        dbc.Nav(
            [
                # HOME
                dbc.NavLink(
                    [
                        html.I(className="fas fa-home"),
                        html.Span("Home", className="navbar_span"),
                    ],
                    href="/",
                    active="exact",
                    className="navbar_entity",
                ),
                # CREATE STUDY
                dbc.NavLink(
                    [
                        html.I(className="fas fa-pen-to-square"),
                        html.Span("Create Study", className="navbar_span"),
                    ],
                    href="/create_study",
                    active="exact",
                    className="navbar_entity",
                ),
                # STUDY DESCRIPTION
                dbc.NavLink(
                    [
                        html.I(className="fas fa-layer-group"),
                        html.Span(
                            "Study Description",
                            className="navbar_span",
                        ),
                    ],
                    href="/study_description",
                    active="exact",
                    className="navbar_entity",
                ),
                # PATWAY
                dbc.NavLink(
                    [
                        html.I(className="fas fa-diagram-project"),
                        html.Span(
                            "Pathways Analysis",
                            className="navbar_span",
                        ),
                    ],
                    href="/pathway_analysis",
                    active="exact",
                    className="navbar_entity",
                ),
                # CLINICAL
                dbc.NavLink(
                    [
                        html.I(className="fas fa-table"),
                        html.Span("Clinical Data", className="navbar_span"),
                    ],
                    href="/clinical_data",
                    active="exact",
                    className="navbar_entity",
                ),
                # SURVIVAL ANALYSIS
                dbc.NavLink(
                    [
                        html.I(className="fa-solid fa-chart-line"),
                        html.Span(
                            "Survival Analysis",
                            className="navbar_span",
                        ),
                    ],
                    href="/survival_analysis",
                    active="exact",
                    className="navbar_entity",
                ),
                # CLUSTER COMP
                dbc.NavLink(
                    [
                        html.I(className="fas fa-code-compare"),
                        html.Span(
                            "Cluster Comparision",
                            className="navbar_span",
                        ),
                    ],
                    href="/cluster_comparision",
                    active="exact",
                    className="navbar_entity",
                ),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    className="sidebar",
)
CONTENT = html.Div(id="page-content", className="content")
APP.layout = html.Div([dcc.Location(id="url", refresh=True), SIDEBAR, CONTENT])

# HOMEPAGE
PAGE_HOME = [
    dbc.Row([html.Img(src=APP_LOGO, className="homepage_logo")]),
    dbc.Row(
        [html.H1(APP_NAME, style={"text-align": "center"})],
    ),
    dbc.Row(
        [
            html.H4(
                [
                    html.Span("ne"),
                    html.Span("T", style={"color": "green"}),
                    html.Span("w"),
                    html.Span("OR", style={"color": "green"}),
                    html.Span("k "),
                    html.Span("T", style={"color": "green"}),
                    html.Span("ool f"),
                    html.Span("O", style={"color": "green"}),
                    html.Span("r mutat"),
                    html.Span("I", style={"color": "green"}),
                    html.Span("onal clu"),
                    html.Span("S", style={"color": "green"}),
                    html.Span("t"),
                    html.Span("E", style={"color": "green"}),
                    html.Span("ring"),
                ],
                style={"text-align": "center"},
            ),
        ],
    ),
    html.Hr(),
    dbc.Row(
        [
            dbc.Col([html.H5("Selected Study")], lg=2),
            dbc.Col([dcc.Dropdown(id="dd-study")], lg=4),
        ],
        justify="center",
    ),
    html.Hr(),
    # HEADER
    dbc.Row(
        [
            html.Span(
                id="dd-study-header",
                className="dd_study_header",
            ),
        ],
    ),
    # PATIENT
    dbc.Row(
        [
            html.Span(
                id="dd-study-info-patient",
                className="dd_study_infos",
            ),
        ],
    ),
    # GENE/VARIANT
    # CLUSTER/MODULARITY
    dbc.Row(
        [
            dbc.Col(
                html.Span(
                    id="dd-study-info-gene",
                    className="dd_study_infos",
                ),
                lg=6,
            ),
            dbc.Col(
                html.Span(
                    id="dd-study-info-variant",
                    className="dd_study_infos",
                ),
                lg=6,
            ),
        ],
    ),
    # CLUSTER/MODULARITY
    dbc.Row(
        [
            dbc.Col(
                html.Span(
                    id="dd-study-info-cluster",
                    className="dd_study_infos",
                ),
                lg=6,
            ),
            dbc.Col(
                html.Span(
                    id="dd-study-info-modularity",
                    className="dd_study_infos",
                ),
                lg=6,
            ),
        ],
    ),
]


# PAGING
@callback(Output("page-content", "children"), Input("url", "pathname"))
def redirect_pages(pathname):
    # If the user tries to reach a different page, return a 404 message
    selected_page = html.Div(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ],
        className="p-3 bg-light rounded-3",
    )
    if pathname == "/":
        selected_page = PAGE_HOME
    elif pathname == "/create_study":
        selected_page = PAGE_CREATE_STUDY
    elif pathname == "/study_description":
        selected_page = PAGE_STUDY_DESCRIPTION
    elif pathname == "/pathway_analysis":
        selected_page = PAGE_PATHWAY_ANALYSIS
    elif pathname == "/cluster_comparision":
        selected_page = PAGE_CLUSTER_COMPARISION
    elif pathname == "/clinical_data":
        selected_page = PAGE_CLINICAL_DATA
    elif pathname == "/survival_analysis":
        selected_page = PAGE_SURVIVAL_ANALYSIS
    return selected_page


@callback(
    Output("dd-study", "options"),
    Output("dd-study", "value"),
    Input("dd-study", "n_clicks"),
)
def update_dropdown_liststudy(_):
    global CONTEXT_DATA
    CONTEXT_DATA["list_studies"] = (
        [] if not Path("study").exists() else os.listdir("study")
    )
    return CONTEXT_DATA["list_studies"], CONTEXT_DATA["name_study"]


# SELECT STUDY
@callback(
    Output("dd-study-header", "children"),
    Output("dd-study-info-patient", "children"),
    Output("dd-study-info-gene", "children"),
    Output("dd-study-info-variant", "children"),
    Output("dd-study-info-cluster", "children"),
    Output("dd-study-info-modularity", "children"),
    Input("dd-study", "value"),
)
def select_study(value) -> str:
    global CONTEXT_DATA
    global GRAPH
    global CLUSTERS_INDEX
    global DF_CLINICAL_DATA
    global NUMERIC_COLUMNS_CLINICAL
    global ALL_COLUMNS_CLINICAL
    global CLUSTER_SELECTED
    global CLUSTER_SELECTED_MULTI
    global BOX_FIG_SELECTED_1
    global BOX_FIG_SELECTED_2

    if value is None:
        return no_update, no_update, no_update, no_update, no_update, no_update
    # prevent reload
    if value == CONTEXT_DATA["name_study"]:
        return (
            f"Study {value}",
            f"Total patient: {CONTEXT_DATA['stats']['num_patient']}",
            f"Total gene: {CONTEXT_DATA['stats']['num_gene']}",
            f"Total variant: {CONTEXT_DATA['stats']['num_variant']}",
            f"Total cluster: {CONTEXT_DATA['stats']['num_cluster']}",
            f"Cluster modularity: {CONTEXT_DATA['stats']['modularity']}",
        )

    # Load config
    with Path("study", value, "config.json").open("r") as json_data:
        d = json.load(json_data)
        CONTEXT_DATA["config"] = d

    CONTEXT_DATA["name_study"] = value
    CONTEXT_DATA["out_root_path"] = Path(
        "study",
        CONTEXT_DATA["name_study"],
        "output",
    )
    GRAPH = np.load(
        CONTEXT_DATA["out_root_path"].joinpath("graph.npy"),
        allow_pickle="TRUE",
    ).item()

    CLUSTERS_INDEX = [int(c) for c in set(GRAPH.vs["cluster"])]

    DF_CLINICAL_DATA = pd.read_csv(
        CONTEXT_DATA["out_root_path"].joinpath("cluster_clinical_data.csv"),
        sep="\t",
        engine="python",
    )
    if len(DF_CLINICAL_DATA.columns) > 1:
        df_clinical_data_all = DF_CLINICAL_DATA.copy()
        df_clinical_data_all["cluster"] = "ALL"
        DF_CLINICAL_DATA = pd.concat([DF_CLINICAL_DATA, df_clinical_data_all])
        NUMERIC_COLUMNS_CLINICAL = list(
            DF_CLINICAL_DATA.select_dtypes(include=np.number).columns,
        )
        ALL_COLUMNS_CLINICAL = list(DF_CLINICAL_DATA.columns)
        ALL_COLUMNS_CLINICAL.remove("cluster")
        BOX_FIG_SELECTED_1 = ALL_COLUMNS_CLINICAL[0]
        BOX_FIG_SELECTED_2 = ALL_COLUMNS_CLINICAL[-1]
    # Reset selection
    CLUSTER_SELECTED = CLUSTERS_INDEX[0]
    CLUSTER_SELECTED_MULTI = ["ALL", CLUSTERS_INDEX[0]]

    # STATS
    stats = pd.read_csv(
        CONTEXT_DATA["out_root_path"].joinpath("numerosity_cluster.csv"),
        sep="\t",
        engine="python",
    )
    CONTEXT_DATA["stats"]["num_patient"] = stats["Patient"].sum()
    CONTEXT_DATA["stats"]["num_gene"] = stats["Gene"].sum()
    CONTEXT_DATA["stats"]["num_variant"] = stats["Variant"].sum()
    CONTEXT_DATA["stats"]["num_cluster"] = len(stats)
    with (
        CONTEXT_DATA["out_root_path"]
        .joinpath("modularity.info")
        .open("r") as f
    ):
        CONTEXT_DATA["stats"]["modularity"] = f.readline()

    return (
        f"Study {value}",
        f"Total patient: {CONTEXT_DATA['stats']['num_patient']}",
        f"Total gene: {CONTEXT_DATA['stats']['num_gene']}",
        f"Total variant: {CONTEXT_DATA['stats']['num_variant']}",
        f"Total cluster: {CONTEXT_DATA['stats']['num_cluster']}",
        f"cluster modularity: {CONTEXT_DATA['stats']['modularity']}",
    )


# CREATE STUDY
PAGE_CREATE_STUDY = [
    dcc.ConfirmDialog(id="confirm-study", message="..."),
    dbc.Row([html.H2("Study Configuration")]),
    # NAME STUDY
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Name Study:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Input(
                        id="input_namestudy",
                        type="text",
                        placeholder="Name study",
                        style={"width": "100%"},
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # NAME STUDY
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Seed trials:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        [5_000, 10_000, 25_000, 50_000, 100_000, 500_000],
                        placeholder=10_000,
                        id="seed-trials",
                        style={"width": "100%"},
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # DATA MUTATIONAL
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("File Data Mutational:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    du.Upload(
                        id="mutational-file",
                        text="Upload Mutational File",
                        chunk_size=100,
                    ),
                ],
                width=3,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        ["\\t", ",", ";"],
                        placeholder="Select separator",
                        id="mutational-separator",
                    ),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Input(
                        id="mutational-skiprow",
                        type="number",
                        placeholder="Skiprows",
                        min=0,
                    ),
                ],
                width=2,
            ),
        ],
    ),
    # COLUMN SAMPLE NAME
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Column Sample Name:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-column-sample-name-mutation",
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # COLUMN GENE NAME
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Column Gene:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-column-gene",
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # IDENTIFIER COLUMNS
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Identifier Columns:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-identifier-columns",
                        multi=True,
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # COLUMN VAF
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Column VAF:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-vaf",
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # COLUMN VAF SCORE
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("VAF Score:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Input(
                        id="input-vaf-score",
                        type="number",
                        placeholder="VAF Score",
                        min=0,
                        max=1,
                        step=0.01,
                        style={"width": "100%"},
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # CLINICAL PATIENT
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("File Clinical Patient:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    du.Upload(
                        id="clinical-patient-file",
                        text="Upload Clinical Patient File",
                        chunk_size=100,
                    ),
                ],
                width=3,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        ["\\t", ",", ";"],
                        placeholder="Select separator",
                        id="clinical-patient-separator",
                    ),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Input(
                        id="clinical-patient-skiprow",
                        type="number",
                        placeholder="Skiprows",
                        min=0,
                    ),
                ],
                width=2,
            ),
        ],
    ),
    # COLUMN PATIENT NAME
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Column Patient Name:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-column-patient-name",
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # COLUMN SURVIVAL EVENT
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Column Survival Event:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-column-survival-event",
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # COLUMN SURVIVAL TIME
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Column Survival Time:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-column-survival-time",
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # CLINICAL SAMPLE
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("File Clinical Sample:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    du.Upload(
                        id="clinical-sample-file",
                        text="Upload Clinical Sample File",
                        chunk_size=100,
                    ),
                ],
                width=3,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        ["\\t", ",", ";"],
                        placeholder="Select separator",
                        id="clinical-sample-separator",
                    ),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Input(
                        id="clinical-sample-skiprow",
                        type="number",
                        placeholder="Skiprows",
                        min=0,
                    ),
                ],
                width=2,
            ),
        ],
    ),
    # COLUMN SAMPLE NAME
    dbc.Row(
        [
            dbc.Col(
                [
                    html.H5("Column Sample Name:"),
                ],
                width=2,
            ),
            dbc.Col(
                [
                    dcc.Dropdown(
                        options=[],
                        id="dd-column-sample-name",
                    ),
                ],
                width=3,
            ),
        ],
    ),
    # CREATE STUDY BUTTON
    dbc.Row(
        [
            html.Button("Create Study", id="create-study-button", n_clicks=0),
        ],
    ),
    # LOADING
    dbc.Row(
        [
            dcc.Loading(
                children=[
                    html.Span(
                        id="create-study-log",
                        className="create_study_log",
                    ),
                ],
            ),
        ],
    ),
]


# dropdown column patient name
@callback(
    Output("dd-column-patient-name", "options"),
    Output("dd-column-survival-event", "options"),
    Output("dd-column-survival-time", "options"),
    [
        Input("clinical-patient-file", "isCompleted"),
        State("clinical-patient-file", "fileNames"),
        Input("clinical-patient-separator", "value"),
        Input("clinical-patient-skiprow", "value"),
    ],
    prevent_initial_call=True,
)
def update_list_columns_patient_name(loaded, filename, sep, skip):
    if not loaded or sep is None or skip is None:
        return [], [], []
    df_patient = pd.read_csv(
        Path("temp", filename[0]),
        sep=sep,
        skiprows=skip,
        engine="python",
        nrows=0,
    )
    return df_patient.columns, df_patient.columns, df_patient.columns


# dropdown column sample name
@callback(
    Output("dd-column-sample-name", "options"),
    [
        Input("clinical-sample-file", "isCompleted"),
        State("clinical-sample-file", "fileNames"),
        Input("clinical-sample-separator", "value"),
        Input("clinical-sample-skiprow", "value"),
    ],
    prevent_initial_call=True,
)
def update_list_columns_sample_name(loaded, filename, sep, skip):
    if not loaded or sep is None or skip is None:
        return []
    df_sample = pd.read_csv(
        Path("temp", filename[0]),
        sep=sep,
        skiprows=skip,
        engine="python",
        nrows=0,
    )
    return df_sample.columns


# dropdown mutation
@callback(
    Output("dd-column-gene", "options"),
    Output("dd-identifier-columns", "options"),
    Output("dd-vaf", "options"),
    Output("dd-column-sample-name-mutation", "options"),
    [
        Input("mutational-file", "isCompleted"),
        State("mutational-file", "fileNames"),
        Input("mutational-separator", "value"),
        Input("mutational-skiprow", "value"),
    ],
    prevent_initial_call=True,
)
def update_list_columns_mutation(loaded, filename, sep, skip):
    if not loaded or sep is None or skip is None:
        return [], [], [], []
    df_mut = pd.read_csv(
        Path("temp", filename[0]),
        sep=sep,
        skiprows=skip,
        engine="python",
        nrows=0,
    )
    return df_mut.columns, df_mut.columns, df_mut.columns, df_mut.columns


@callback(
    Output("confirm-study", "displayed"),
    Output("confirm-study", "message"),
    Input("create-study-button", "n_clicks"),
    State("input_namestudy", "value"),
    State("mutational-file", "fileNames"),
    State("mutational-separator", "value"),
    State("mutational-skiprow", "value"),
    State("clinical-patient-file", "fileNames"),
    State("clinical-patient-separator", "value"),
    State("clinical-patient-skiprow", "value"),
    State("clinical-sample-file", "fileNames"),
    State("clinical-sample-separator", "value"),
    State("clinical-sample-skiprow", "value"),
    State("dd-column-patient-name", "value"),
    State("dd-column-survival-event", "value"),
    State("dd-column-survival-time", "value"),
    State("dd-column-sample-name", "value"),
    State("dd-column-sample-name-mutation", "value"),
    State("dd-column-gene", "value"),
    State("dd-identifier-columns", "value"),
    State("dd-vaf", "value"),
    State("input-vaf-score", "value"),
    State("seed-trials", "value"),
    prevent_initial_call=True,
)
def create_study(
    _,
    input_namestudy,
    mutational_filename,
    mutational_separator,
    mutational_skiprow,
    clinical_patient_filename,
    clinical_patient_separator,
    clinical_patient_skiprow,
    clinical_sample_filename,
    clinical_sample_separator,
    clinical_sample_skiprow,
    c_patient_name,
    c_surv_event,
    c_surv_time,
    c_sample_name,
    c_sample_mutation,
    c_gene,
    c_identifier,
    c_vaf,
    vaf_score,
    seed_trials
):
    global CONTEXT_DATA
    global PATH_CONFIG
    global READY_FOR_CREATION

    READY_FOR_CREATION = False
    # CHECK INPUT DATA
    if input_namestudy is None:
        return True, "Please insert name for study"
    CONTEXT_DATA["name_study_input"] = input_namestudy
    if Path("study", input_namestudy, "input").exists():
        return True, "Name study already exist"
    if mutational_filename is None:
        return True, "Please insert file for mutational data"
    if mutational_separator is None:
        return True, "Please select separator for mutational data"
    if mutational_skiprow is None:
        return True, "Please select skiprow for mutational data"
    if c_gene is None:
        return True, "Please select column for gene"
    if c_sample_mutation is None:
        return (
            True,
            "Please select column for sample name on mutation file",
        )
    if c_identifier == []:
        return (
            True,
            "Please select at least one column for mutation identifier",
        )
    if clinical_patient_filename is not None:
        if clinical_patient_separator is None:
            return True, "Please select separator for clinical patient data"
        if clinical_patient_skiprow is None:
            return True, "Please select skiprow for clinical patient data"
        if c_patient_name is None:
            return True, "Please select column for patient name"
        if (c_surv_event is None and c_surv_time is not None) or (
            c_surv_event is not None and c_surv_time is None
        ):
            return True, "Please select survival envent and survival time"
    if clinical_sample_filename is not None:
        if clinical_sample_separator is None:
            return True, "Please select separator for clinical sample data"
        if clinical_sample_skiprow is None:
            return True, "Please select skiprow for clinical sample data"
        if c_sample_name is None:
            return (
                True,
                "Please select column for sample name on clinical sample file",
            )
    if seed_trials is None:
        seed_trials = 10_000
    # CREATE STUDY FOLDER
    Path("study", input_namestudy, "input").mkdir(parents=True, exist_ok=True)
    # MUTATIONAL DATA
    Path("temp", mutational_filename[0]).rename(
        Path("study", input_namestudy, "input", "mutational_data.txt"),
    )
    dialog_message = f"Create study '{input_namestudy}'?"
    # CLINICAL PATIENT
    if clinical_patient_filename is not None:
        Path("temp", clinical_patient_filename[0]).rename(
            Path("study", input_namestudy, "input", "clinical_patient.txt"),
        )
    # CLINICAL SAMPLE
    if clinical_sample_filename is not None:
        Path("temp", clinical_sample_filename[0]).rename(
            Path("study", input_namestudy, "input", "clinical_sample.txt"),
        )
    # GENERATE JSON CONFIG
    config_dict = {}
    config_dict["paths"] = {}
    config_dict["clinical_data"] = {}
    config_dict["mutation"] = {}
    config_dict["name"] = input_namestudy
    config_dict["paths"]["data_mutational"] = str(
        Path(
            "study",
            input_namestudy,
            "input",
            "mutational_data.txt",
        ),
    )
    config_dict["paths"]["data_mutational_sep"] = "\t"
    config_dict["paths"]["data_mutational_skip"] = mutational_skiprow
    config_dict["paths"]["data_clinical_patient"] = ""
    config_dict["paths"]["data_clinical_sample_sep"] = "\t"
    config_dict["paths"]["data_clinical_sample_skip"] = 0
    config_dict["paths"]["data_clinical_sample"] = ""
    config_dict["paths"]["data_clinical_patient_sep"] = "\t"
    config_dict["paths"]["data_clinical_patient_skip"] = 0
    config_dict["mutation"]["column_gene"] = c_gene
    config_dict["mutation"]["column_sample_name"] = c_sample_mutation
    config_dict["seed_trials"] = seed_trials
    # REMOVE EXTRA SEPARATOR BEFOR JOIN
    c_identifier_clean = [col.replace(";", "") for col in c_identifier]
    config_dict["mutation"]["identifier_columns"] = ";".join(
        c_identifier_clean,
    )
    config_dict["clinical_data"]["column_patient_name"] = ""
    config_dict["clinical_data"]["column_sample_name"] = ""
    config_dict["clinical_data"]["column_surv_event"] = ""
    config_dict["clinical_data"]["column_surv_time"] = ""
    config_dict["mutation"]["vaf_score"] = vaf_score
    config_dict["mutation"]["vaf_column"] = "" if c_vaf is None else c_vaf
    if clinical_patient_filename is not None:
        config_dict["paths"]["data_clinical_patient"] = str(
            Path(
                "study",
                input_namestudy,
                "input",
                "clinical_patient.txt",
            ),
        )
        config_dict["paths"]["data_clinical_patient_skip"] = (
            clinical_patient_skiprow
        )
        config_dict["clinical_data"]["column_patient_name"] = c_patient_name
        if c_surv_time and c_surv_event:
            config_dict["clinical_data"]["column_surv_event"] = c_surv_event
            config_dict["clinical_data"]["column_surv_time"] = c_surv_time
    if clinical_sample_filename is not None:
        config_dict["paths"]["data_clinical_sample"] = str(
            Path(
                "study",
                input_namestudy,
                "input",
                "clinical_sample.txt",
            ),
        )
        config_dict["paths"]["data_clinical_sample_skip"] = (
            clinical_sample_skiprow
        )
        config_dict["clinical_data"]["column_sample_name"] = c_sample_name
    PATH_CONFIG = Path("study", input_namestudy, "config.json")
    with PATH_CONFIG.open("w", encoding="utf-8") as f:
        json.dump(config_dict, f, indent=4)

    READY_FOR_CREATION = True
    return True, dialog_message


@callback(
    Output("create-study-log", "children"),
    Output("confirm-study", "submit_n_clicks"),
    Output("confirm-study", "cancel_n_clicks"),
    [
        Input("confirm-study", "submit_n_clicks"),
        Input("confirm-study", "cancel_n_clicks"),
    ],
    prevent_initial_call=True,
)
def confirm_study_prompt(submit_n_clicks, _):
    global CONTEXT_DATA
    global PATH_CONFIG

    if not READY_FOR_CREATION:
        return html.Div(""), 0, 0

    if submit_n_clicks:
        with conversion.localconverter(default_converter):
            os.popen(f"python tortoise.py -c {PATH_CONFIG}").read()
            select_study(CONTEXT_DATA["name_study_input"])
        return (
            html.Div(f"Study '{CONTEXT_DATA['name_study_input']}' created!"),
            0,
            0,
        )
    # Delete temp files
    study_path = Path("study", CONTEXT_DATA["name_study_input"])
    for x in study_path.glob("*"):
        if x.is_dir():
            for y in x.glob("*"):
                y.unlink()
            x.rmdir()
        else:
            x.unlink()
    study_path.rmdir()
    PATH_CONFIG = None
    CONTEXT_DATA["name_study_input"] = None

    return html.Div(""), 0, 0


# STUDY DESCRIPTION
PAGE_STUDY_DESCRIPTION = [
    # DIALOG INFO NODE
    dbc.Modal(
        [],
        id="modal-lg",
        size="lg",
        is_open=False,
    ),
    # ROW LINE PARAMS
    dbc.Row(
        [
            # CLUSTER SELECTOR
            dbc.Col(
                [
                    html.Span("Cluster selected", className="span_selector"),
                    dcc.Dropdown(id="dd-cluster"),
                    html.Br(),
                ],
                lg=6,
            ),
            html.Hr(),
        ],
    ),
    # 1 ROW
    dbc.Row(
        [
            dbc.Col(
                [
                    html.Span("Number of patients"),
                    html.Hr(),
                    html.Span(0, id="span_n_patient"),
                ],
                width=2,
                className="info_block add-border",
            ),
            dbc.Col(
                [
                    html.Span("Number of variants"),
                    html.Hr(),
                    html.Span(0, id="span_n_variants"),
                ],
                width=2,
                className="info_block add-border",
            ),
            dbc.Col(
                [
                    html.Span("Number of genes"),
                    html.Hr(),
                    html.Span(0, id="span_n_genes"),
                ],
                width=2,
                className="info_block add-border",
            ),
            dbc.Col(
                [
                    html.Span("Variant centroid"),
                    html.Hr(),
                    html.Span("None", id="span_variant_centroids"),
                ],
                width=2,
                className="info_block add-border",
            ),
        ],
        justify="evenly",
    ),
    dbc.Row([html.Br()]),
    # 2 ROW
    dbc.Row(
        [
            # CYTOSCAPE BLOCK
            dbc.Col(
                [
                    # CYTOSCAPE GRAPH
                    cyto.Cytoscape(
                        className="add-border",
                        id="cytoscape-graph",
                        style={"width": "100%", "height": "35vh"},
                        stylesheet=[
                            # NOME SOPRA
                            {
                                "selector": "node",
                                "style": {
                                    "content": "data(name)",
                                    "font-size": "5px",
                                },
                            },
                            # PAZIENTI TRIANGOLI ROSSI
                            {
                                "selector": ".PATIENT",
                                "style": {
                                    "background-color": "coral",
                                    "shape": "triangle",
                                },
                            },
                            # VARIANTI CERCHI BLUE
                            {
                                "selector": ".VARIANT",
                                "style": {
                                    "background-color": "royalblue",
                                    "shape": "circle",
                                },
                            },
                            # SELECTEDƒco
                            {
                                "selector": ":selected",
                                "style": {
                                    "background-color": "#02cd79",
                                },
                            },
                        ],
                        minZoom=0.1,
                        maxZoom=2,
                        responsive=True,
                    ),
                    # LAYOUT SELECTOR
                    dcc.RadioItems(
                        options=[
                            "cose",
                            "concentric",
                            "grid",
                            "circle",
                            "breadthfirst",
                            "klay",
                        ],
                        value="cose",
                        inline=True,
                        id="radio-layouts",
                        persistence=True,
                        persistence_type="memory",
                        style={"width": "100%", "height": "2vh"},
                    ),
                ],
                lg=6,
            ),
            # FIGURE PIE
            dbc.Col(
                [
                    dcc.Graph(
                        id="fig_pie",
                        className="add-border",
                        style={"width": "100%", "height": "35vh"},
                    ),
                ],
                lg=6,
            ),
        ],
        justify="evenly",
    ),
    dbc.Row([html.Br()]),
    # 3 ROW
    dbc.Row(
        [
            # FIGURE DEGREE
            dbc.Col(
                [
                    dcc.Graph(
                        id="fig_degree",
                        className="add-border",
                        style={"width": "100%", "height": "35vh"},
                    ),
                ],
                lg=8,
            ),
        ],
        justify="evenly",
    ),
]


@callback(
    Output("dd-cluster", "options"),
    Output("dd-cluster", "value"),
    Input("dd-cluster", "n_clicks"),
)
def dropdpwn_cluster(_):
    return CLUSTERS_INDEX, CLUSTER_SELECTED


# SELECT SINGLE NODE GRAPH
@callback(
    Output("modal-lg", "children"),
    Output("modal-lg", "is_open"),
    Input("cytoscape-graph", "tapNodeData"),
    prevent_initial_call=True,
)
def display_node_data(data_dict):
    temp = ""
    if data_dict["vertex_type"] == "VARIANT":
        term_included = ["name", "gene", "sost_amm"]
        temp = ""
        for k, v in data_dict.items():
            if k in term_included:
                temp += f"**{k}**:{v}\n"
    else:
        term_excluded = [
            "vertex_type",
            "variants",
            "color_vertex",
            "shape_vertex",
            "gene",
            "sost_amm",
            "consequence",
            "color",
            "cluster",
            "id",
            "timeStamp",
        ]
        temp = ""
        for k, v in data_dict.items():
            if k not in term_excluded:
                temp += f"**{k}**:{v}\n"
    return [
        dbc.ModalHeader(dbc.ModalTitle(data_dict["name"])),
        dcc.Markdown(temp, className="markdown"),
    ], True


# UPDATE CLUSTER LAYOUT
@callback(Output("cytoscape-graph", "layout"), Input("radio-layouts", "value"))
def update_graph(layout):
    return {"name": layout}


# SELECT CLUSTER INDEX
@callback(
    Output("cytoscape-graph", "elements"),
    Output("fig_pie", "figure"),
    Output("fig_degree", "figure"),
    Output("span_n_patient", "children"),
    Output("span_n_variants", "children"),
    Output("span_n_genes", "children"),
    Output("span_variant_centroids", "children"),
    Input("dd-cluster", "value"),
)
def update_cluster(cluster):
    if cluster is None:
        return (
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
        )
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    # CLUSTER ELEMENTS
    cluster_elements = filter_graph(cluster)
    # FIGURE PIE
    df_gene = pd.read_csv(
        Path(
            CONTEXT_DATA["out_root_path"],
            "Gene_Count",
            f"genes_cluster_{cluster}.csv",
        ),
        sep="\t",
        engine="python",
    )
    fig_pie = px.pie(
        df_gene,
        values="COUNT",
        names="GENE",
        title="Number Mutation for Gene",
    )
    fig_pie.update_traces(textposition="inside", textinfo="label")
    # VARIANT NUMBERS
    df_variant = pd.read_csv(
        Path(
            CONTEXT_DATA["out_root_path"],
            "Variants_Degree",
            f"variants_degree_cluster{cluster}.csv",
        ),
        sep="\t",
        engine="python",
    )
    n_variants = len(df_variant)
    # FIGURE DEGREE
    df_variant = df_variant.sort_values(by=["Degree"], ascending=False)[:15]
    fig_degree = px.bar(
        df_variant,
        x="Variants",
        y="Degree",
        title="Mutation Degree",
    )
    # PATIENTS NUMBER
    n_patients = len(
        [
            1
            for e in cluster_elements
            if e["data"].get("vertex_type", "") == "PATIENT"
        ],
    )
    # GENE NUMBERS
    n_genes = len(df_gene)
    # VARIANT CENTROID
    if (
        n_variants > 1
        and df_variant.iloc[0]["Degree"] == df_variant.iloc[1]["Degree"]
    ):
        variant_centroids = "More than one"
    else:
        variant_centroids = df_variant.iloc[0]["Variants"]
    # RETURN
    return (
        cluster_elements,
        fig_pie,
        fig_degree,
        n_patients,
        n_variants,
        n_genes,
        variant_centroids,
    )


# PATHWAY ANALYSIS
PAGE_PATHWAY_ANALYSIS = [
    # ROW LINE PARAMS
    dbc.Row(
        [
            # CLUSTER SELECTOR
            dbc.Col(
                [
                    html.Span("Cluster selected", className="span_selector"),
                    dcc.Dropdown(id="dd-cluster"),
                    html.Br(),
                ],
                lg=6,
            ),
            # PVALUE SELECTOR
            dbc.Col(
                [
                    html.Span("PValue threshold", className="span_selector"),
                    dcc.Dropdown(
                        [0.01, 0.05],
                        0.05,
                        id="dd-pvalue",
                        persistence=True,
                        persistence_type="memory",
                    ),
                    html.Br(),
                ],
                lg=3,
            ),
            # ADJUSTED PVALUE SELECTOR
            dbc.Col(
                [
                    html.Span(
                        "Use adjusted PValue",
                        className="span_selector",
                    ),
                    dcc.Dropdown(
                        ["True", "False"],
                        "False",
                        id="dd-adjusted-pvalue",
                        persistence=True,
                        persistence_type="memory",
                    ),
                    html.Br(),
                ],
                lg=3,
            ),
            html.Hr(),
        ],
    ),
    # FIRST ROW
    dbc.Row(
        [
            # FIGURE GO
            dbc.Col(
                [
                    dcc.Graph(
                        id="fig_go",
                        className="add-border",
                        style={"width": "100%", "height": "41vh"},
                    ),
                    dcc.RadioItems(
                        options=[
                            {
                                "label": "Biological Function",
                                "value": "biological",
                            },
                            {
                                "label": "Molecular Function",
                                "value": "molecular",
                            },
                            {
                                "label": "Cellular Component",
                                "value": "cellular",
                            },
                        ],
                        # Valore predefinito
                        value="biological",
                        labelStyle={"display": "inline-block"},
                        id="radio_fig_go",
                        persistence=True,
                        persistence_type="memory",
                    ),
                ],
                lg=6,
            ),
            # FIGURE KEGG
            dbc.Col(
                [
                    dcc.Graph(
                        id="fig_kegg",
                        className="add-border",
                        style={"width": "100%", "height": "41vh"},
                    ),
                ],
                lg=6,
            ),
            # FIGURE REACTOME
            dbc.Col(
                [
                    dcc.Graph(
                        id="fig_reactome",
                        className="add-border",
                        style={"width": "100%", "height": "41vh"},
                    ),
                ],
                lg=6,
            ),
            # FIGURE WIKI
            dbc.Col(
                [
                    dcc.Graph(
                        id="fig_wiki",
                        className="add-border",
                        style={"width": "100%", "height": "41vh"},
                    ),
                ],
                lg=6,
            ),
        ],
    ),
]


def generate_pathway_fig(df, pvalue, adj_pvalue, title):
    col_name = "P.value"
    label_name = "Pvalue"
    if adj_pvalue == "True":
        col_name = "Adjusted.P.value"
        label_name = "Adjusted Pvalue"
    df_data = df[df[col_name] < pvalue]
    df_data = df_data.sort_values(by=[col_name], ascending=False)[-25:]
    df_data["Count_gene"] = df_data.apply(
        lambda e: len(e["Genes"].split(";")),
        axis=1,
    )
    fig = px.bar(
        df_data,
        x="Count_gene",
        y="Term",
        hover_data=["Overlap"],
        color=col_name,
        title=title,
        labels={col_name: label_name},
    )
    fig.update_layout(
        xaxis_title="Genes_Count",
        yaxis_title="Terms",
        legend_title=label_name,
        coloraxis_colorbar={"exponentformat": "e"},
    )
    return fig


# UPDATE GO FIGURE
@callback(
    Output("fig_go", "figure"),
    [
        Input("dd-cluster", "value"),
        Input("dd-pvalue", "value"),
        Input("dd-adjusted-pvalue", "value"),
        Input("radio_fig_go", "value"),
    ],
)
def update_go(cluster, pvalue, adj_pvalue, p_type):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df_data = pd.read_csv(
        Path(
            CONTEXT_DATA["out_root_path"],
            "Arricchimento_all_genes",
            "GO",
            f"{p_type}_{cluster}.csv",
        ),
        engine="python",
    )
    return generate_pathway_fig(df_data, pvalue, adj_pvalue, "GO")


# UPDATE KEGG FIGURE
@callback(
    Output("fig_kegg", "figure"),
    [
        Input("dd-cluster", "value"),
        Input("dd-pvalue", "value"),
        Input("dd-adjusted-pvalue", "value"),
    ],
)
def update_kegg(cluster, pvalue, adj_pvalue):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df_data = pd.read_csv(
        Path(
            CONTEXT_DATA["out_root_path"],
            "Arricchimento_all_genes",
            "KEGG",
            f"kegg_{cluster}.csv",
        ),
        engine="python",
    )
    return generate_pathway_fig(df_data, pvalue, adj_pvalue, "KEGG")


# UPDATE REACTOME FIGURE
@callback(
    Output("fig_reactome", "figure"),
    [
        Input("dd-cluster", "value"),
        Input("dd-pvalue", "value"),
        Input("dd-adjusted-pvalue", "value"),
    ],
)
def update_reactome(cluster, pvalue, adj_pvalue):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df_data = pd.read_csv(
        Path(
            CONTEXT_DATA["out_root_path"],
            "Arricchimento_all_genes",
            "REACTOME",
            f"reactome_{cluster}.csv",
        ),
        engine="python",
    )
    return generate_pathway_fig(df_data, pvalue, adj_pvalue, "REACTOME")


# UPDATE WIKI FIGURE
@callback(
    Output("fig_wiki", "figure"),
    [
        Input("dd-cluster", "value"),
        Input("dd-pvalue", "value"),
        Input("dd-adjusted-pvalue", "value"),
    ],
)
def update_wiki(cluster, pvalue, adj_pvalue):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df_data = pd.read_csv(
        Path(
            CONTEXT_DATA["out_root_path"],
            "Arricchimento_all_genes",
            "WIKI",
            f"wiki_{cluster}.csv",
        ),
        engine="python",
    )
    return generate_pathway_fig(df_data, pvalue, adj_pvalue, "WIKI")


# CLINICAL DATA
PAGE_CLINICAL_DATA = [
    # DIALOG INFO NODE
    dbc.Modal(
        [],
        id="modal-lg",
        size="lg",
        is_open=False,
    ),
    # ROW LINE PARAMS
    dbc.Row(
        [
            # CLUSTER SELECTOR
            dbc.Col(
                [
                    html.Span("Cluster selected", className="span_selector"),
                    dcc.Dropdown(id="dd-cluster"),
                    html.Br(),
                ],
                lg=6,
            ),
            html.Hr(),
        ],
    ),
    dbc.Row(
        [
            # FIGURE BOX_PLOT_1
            dbc.Col(
                [
                    dcc.Dropdown(id="dd-box-1"),
                    dcc.Graph(
                        id="fig_box_plot_1",
                        className="add-border",
                        style={"width": "100%", "height": "40vh"},
                    ),
                ],
                lg=6,
            ),
            # FIGURE BOX_PLOT_2
            dbc.Col(
                [
                    dcc.Dropdown(id="dd-box-2"),
                    dcc.Graph(
                        id="fig_box_plot_2",
                        className="add-border",
                        style={"width": "100%", "height": "40vh"},
                    ),
                ],
                lg=6,
            ),
        ],
    ),
    dbc.Row([dash_table.DataTable(id="table_clinical_data")]),
]


@callback(
    Output("dd-box-1", "options"),
    Output("dd-box-1", "value"),
    Input("dd-box-1", "n_clicks"),
)
def dropdown_box_fig_1(_):
    return ALL_COLUMNS_CLINICAL, BOX_FIG_SELECTED_1


@callback(
    Output("dd-box-2", "options"),
    Output("dd-box-2", "value"),
    Input("dd-box-2", "n_clicks"),
)
def dropdown_box_fig_2(_):
    return ALL_COLUMNS_CLINICAL, BOX_FIG_SELECTED_2


# SINGLE IMAGE BOX/PIE
def func_single_plot(cluster, col_name):
    cluster_values = DF_CLINICAL_DATA[DF_CLINICAL_DATA["cluster"] == cluster][
        col_name
    ]

    if col_name in NUMERIC_COLUMNS_CLINICAL:
        return px.box(cluster_values, y=col_name)

    _temp_dict = dict(cluster_values.value_counts())
    _temp_df = pd.DataFrame(
        {col_name: _temp_dict.keys(), "count": _temp_dict.values()},
    )
    return px.pie(_temp_df, values="count", names=col_name)


# UPDATE BOX_PLOT_1
@callback(
    Output("fig_box_plot_1", "figure"),
    [Input("dd-cluster", "value"), Input("dd-box-1", "value")],
)
def update_box_1(cluster, col_name):
    if cluster is None or col_name is None:
        return no_update
    global BOX_FIG_SELECTED_1
    global CLUSTER_SELECTED
    BOX_FIG_SELECTED_1 = col_name
    CLUSTER_SELECTED = cluster
    return func_single_plot(cluster, col_name)


# UPDATE BOX_PLOT_2
@callback(
    Output("fig_box_plot_2", "figure"),
    [Input("dd-cluster", "value"), Input("dd-box-2", "value")],
)
def update_box_2(cluster, col_name):
    if cluster is None or col_name is None:
        return no_update
    global BOX_FIG_SELECTED_2
    global CLUSTER_SELECTED
    BOX_FIG_SELECTED_2 = col_name
    CLUSTER_SELECTED = cluster
    return func_single_plot(cluster, col_name)


# TABLE CLINICAL_DATA:
@callback(Output("table_clinical_data", "data"), Input("dd-cluster", "value"))
def update_table_clinical_data(cluster):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    cluster_values = DF_CLINICAL_DATA[DF_CLINICAL_DATA["cluster"] == cluster]
    return cluster_values.to_dict("records")


# SURVIVAL ANALYSIS
PAGE_SURVIVAL_ANALYSIS = [
    # DIALOG INFO NODE
    dbc.Modal(
        [],
        id="modal-lg",
        size="lg",
        is_open=False,
    ),
    # ROW LINE PARAMS
    dbc.Row(
        [
            # CLUSTER SELECTOR
            dbc.Col(
                [
                    html.Span("Cluster selected", className="span_selector"),
                    dcc.Dropdown(id="dd-cluster"),
                    html.Br(),
                ],
                lg=6,
            ),
            html.Hr(),
        ],
    ),
    dbc.Row(
        [
            dbc.Col(
                [
                    dcc.Graph(
                        id="survival_figure",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    ),
                ],
                lg=8,
            ),
            dbc.Col(
                [
                    dcc.Graph(
                        id="survival_figure_stat",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    ),
                ],
                lg=4,
            ),
        ],
    ),
    dbc.Row(
        [
            # CLUSTER SELECTOR MULTI
            dbc.Col(
                [
                    html.Span(
                        "Multi Cluster Selection",
                        className="span_selector",
                    ),
                    dcc.Dropdown(
                        id="dd-cluster-multi",
                        multi=True,
                    ),
                    html.Br(),
                ],
                lg=4,
            ),
            html.Hr(),
        ],
    ),
    dbc.Row(
        [
            dbc.Col(
                [
                    dcc.Graph(
                        id="survival_figure_comparison",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    ),
                ],
                lg=8,
            ),
            dbc.Col(
                [
                    dcc.Graph(
                        id="table_test_survival",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    ),
                ],
                lg=4,
            ),
        ],
    ),
]


@callback(
    Output("dd-cluster-multi", "options"),
    Output("dd-cluster-multi", "value"),
    Input("dd-cluster-multi", "n_clicks"),
)
def dropdpwn_multi_cluster(_):
    return ["ALL", *CLUSTERS_INDEX], CLUSTER_SELECTED_MULTI


# SURVIVAL_PLOT
@callback(
    Output("survival_figure", "figure"),
    Output("survival_figure_stat", "figure"),
    Input("dd-cluster", "value"),
)
def update_overall_survival(cluster):
    global CLUSTER_SELECTED
    if cluster is None:
        return no_update, no_update
    CLUSTER_SELECTED = cluster
    col_surv_time = CONTEXT_DATA["config"]["clinical_data"]["column_surv_time"]
    col_surv_event = CONTEXT_DATA["config"]["clinical_data"][
        "column_surv_event"
    ]
    data = DF_CLINICAL_DATA[DF_CLINICAL_DATA["cluster"] == cluster]
    data = data.dropna(subset=[col_surv_time, col_surv_event])
    if len(data) <= 0:
        return no_update, no_update
    data[col_surv_event] = data[col_surv_event].replace({"Yes": 1, "No": 0})
    kmf = KaplanMeierFitter()
    kmf.fit(
        data[col_surv_time].values,
        event_observed=data[col_surv_event].values,
    )
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=kmf.confidence_interval_.index,
            y=kmf.confidence_interval_["KM_estimate_upper_0.95"],
            mode="lines",
            line={"shape": "hv", "width": 0},
            showlegend=False,
        ),
    )
    fig.add_trace(
        go.Scatter(
            x=kmf.confidence_interval_.index,
            y=kmf.confidence_interval_["KM_estimate_lower_0.95"],
            mode="lines",
            line={"shape": "hv", "width": 0},
            fill="tonexty",
            fillcolor="rgb(153,204,255)",
            showlegend=False,
        ),
    )
    fig.update_layout(
        title=f"Kaplan-Meier Curve Cluster: {cluster}",
        xaxis_title="Duration",
        yaxis_title="Survival probability",
        font_size=14,
        xaxis_title_font_size=18,
        yaxis_title_font_size=18,
    )
    fig.add_trace(
        go.Scatter(
            x=kmf.survival_function_.index,
            y=kmf.survival_function_["KM_estimate"],
            line={"shape": "hv", "width": 3, "color": "rgb(0,0,128)"},
            mode="lines",
            showlegend=False,
        ),
    )
    # STAT PIE
    fig_stat = func_single_plot(cluster, col_surv_event)
    fig_stat.update_layout(
        {
            "title": f"Vital Status Cluster: {cluster}",
            "legend_title": "Alive",
        },
    )
    return fig, fig_stat


@callback(
    Output("survival_figure_comparison", "figure"),
    Output("table_test_survival", "figure"),
    Input("dd-cluster-multi", "value"),
)
def update_survival_comparison(list_clusters):
    global CLUSTER_SELECTED_MULTI
    min_cluster = 2
    if len(list_clusters) < min_cluster:
        return no_update, no_update
    CLUSTER_SELECTED_MULTI = list_clusters
    fig = go.Figure()
    fig_stats = None
    kmf = KaplanMeierFitter()
    col_surv_time = CONTEXT_DATA["config"]["clinical_data"]["column_surv_time"]
    col_surv_event = CONTEXT_DATA["config"]["clinical_data"][
        "column_surv_event"
    ]
    df_clinical = DF_CLINICAL_DATA[
        DF_CLINICAL_DATA["cluster"].isin(list_clusters)
    ]
    df_clinical = df_clinical.dropna(subset=[col_surv_time, col_surv_event])
    df_clinical[col_surv_event] = df_clinical[col_surv_event].replace(
        {"Yes": 1, "No": 0},
    )

    for cluster in list(df_clinical["cluster"].unique()):
        cluster_data = df_clinical[df_clinical["cluster"] == cluster]
        if len(cluster_data) == 0:
            list_clusters.remove(cluster)
        else:
            kmf.fit(
                cluster_data[col_surv_time],
                event_observed=cluster_data[col_surv_event],
                label=f"Cluster {cluster}",
            )
            fig.update_layout(
                title="Kaplan-Meier Curve",
                xaxis_title="Duration",
                yaxis_title="Survival probability",
                font_size=14,
                xaxis_title_font_size=18,
                yaxis_title_font_size=18,
            )
            fig.add_trace(
                go.Scatter(
                    x=kmf.survival_function_.index,
                    y=kmf.survival_function_[f"Cluster {cluster}"],
                    line={"shape": "hv", "width": 3},
                    mode="lines",
                    name=f"Cluster {cluster}",
                    showlegend=True,
                ),
            )

    if len(list_clusters) <= 1:
        return no_update, no_update

    df_clinical["cluster_str"] = df_clinical["cluster"].astype(str)
    # Test log-rank
    results = pairwise_logrank_test(
        df_clinical[col_surv_time],
        df_clinical["cluster_str"],
        df_clinical[col_surv_event],
    )
    df_res = results.summary
    df_res = df_res.rename(columns={"-log2(p)": "log2_p"})
    df_res = df_res.reset_index()
    df_res["comparison"] = (
        df_res["level_0"].astype(str) + " " + df_res["level_1"].astype(str)
    )
    # remove col level_0 e level_1
    df_res = df_res.drop(columns=["level_0", "level_1"])
    fig_stats = go.Figure(
        data=[
            go.Table(
                header={
                    "values": list(df_res.columns),
                    "fill_color": "paleturquoise",
                    "align": "left",
                },
                cells={
                    "values": [
                        round(df_res.test_statistic, 4),
                        round(df_res.p, 4),
                        round(df_res.log2_p, 4),
                        df_res.comparison,
                    ],
                    "fill_color": "lavender",
                    "align": "left",
                },
            ),
        ],
    )
    return fig, fig_stats


# CLUSTER COMPARISION
PAGE_CLUSTER_COMPARISION = [
    # ROW LINE PARAMS
    dbc.Row(
        [
            # CLUSTER SELECTOR MULTI
            dbc.Col(
                [
                    html.Span("Cluster selected", className="span_selector"),
                    dcc.Dropdown(id="dd-cluster-multi", multi=True),
                    html.Br(),
                ],
                lg=4,
            ),
            html.Hr(),
        ],
    ),
    # COL VENN
    dbc.Row(
        [
            dbc.Col(
                [
                    html.Img(
                        id="plot-venn",
                        className="add-border",
                        style={"width": "100%", "height": "80vh"},
                    ),
                ],
                lg=6,
            ),
            dbc.Col(
                [
                    dcc.Graph(
                        id="tbl_g_common",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    ),
                    dcc.Dropdown(id="dd-box-1"),
                    dcc.Graph(
                        id="fig_multi_fig",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    ),
                ],
                lg=6,
            ),
        ],
    ),
]


# UPDATE COMPARISION
@callback(Output("plot-venn", "src"), Input("dd-cluster-multi", "value"))
def update_venn(list_clusters):
    global CLUSTER_SELECTED_MULTI
    min_cluster = 2
    if len(list_clusters) < min_cluster:
        return no_update
    CLUSTER_SELECTED_MULTI = list_clusters

    cluster_gene_list = []
    for index in list_clusters:
        gene_values = []
        if index == "ALL":
            gene_values = pd.read_csv(
                Path(
                    CONTEXT_DATA["out_root_path"],
                    "distribution_gene_cluster.csv",
                ),
                sep="\t",
                engine="python",
            )["Gene"].unique()
        else:
            gene_values = pd.read_csv(
                Path(
                    CONTEXT_DATA["out_root_path"],
                    "Gene_Count",
                    f"genes_cluster_{index}.csv",
                ),
                sep="\t",
                engine="python",
            )["GENE"].to_numpy()
        cluster_gene_list.append(gene_values)
    labels = venn.get_labels(cluster_gene_list, fill=["number"])

    fig, ax = None, None
    match len(list_clusters):
        case 2:
            fig, ax = venn.venn2(labels, names=list_clusters, figsize=(10, 10))
        case 3:
            fig, ax = venn.venn3(labels, names=list_clusters, figsize=(10, 10))
        case 4:
            fig, ax = venn.venn4(labels, names=list_clusters, figsize=(10, 10))
        case 5:
            fig, ax = venn.venn5(labels, names=list_clusters, figsize=(10, 10))
        case _:
            return no_update

    ax.set_title("Gene Comparision")
    # SAVE TO BUFFER
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    fig_data = base64.b64encode(buf.getbuffer()).decode("utf-8")
    fig_bar_matplotlib = f"data:image/png;base64,{fig_data}"
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()
    return fig_bar_matplotlib


# TABLE GENES_COMMON:
@callback(Output("tbl_g_common", "figure"), Input("dd-cluster-multi", "value"))
def update_genes_common(list_clusters):
    global CLUSTER_SELECTED_MULTI
    min_cluster = 2
    if len(list_clusters) < min_cluster:
        return no_update
    CLUSTER_SELECTED_MULTI = list_clusters

    cluster_gene_list = {}
    all_gene_set = set()
    for index in list_clusters:
        gene_values = []
        if index == "ALL":
            gene_values = list(
                pd.read_csv(
                    Path(
                        CONTEXT_DATA["out_root_path"],
                        "distribution_gene_cluster.csv",
                    ),
                    sep="\t",
                    engine="python",
                )["Gene"].unique(),
            )
            cluster_gene_list["ALL"] = gene_values
        else:
            gene_values = list(
                pd.read_csv(
                    Path(
                        CONTEXT_DATA["out_root_path"],
                        "Gene_Count",
                        f"genes_cluster_{index}.csv",
                    ),
                    sep="\t",
                    engine="python",
                )["GENE"].unique(),
            )
            cluster_gene_list[index] = gene_values
        all_gene_set.update(gene_values)

    df_genes = pd.DataFrame({"gene": list(all_gene_set)})
    for cluster, genes in cluster_gene_list.items():
        df_genes[cluster] = df_genes["gene"].apply(
            lambda x, g=genes: "🟢" if x in g else "🔴",
        )

    df_genes = df_genes.sort_values(by=list_clusters, ascending=False)
    values = [df_genes["gene"].to_numpy()]
    c_name = ["Gene"]
    for i in list_clusters:
        c_name.append(f"Cluster {i}")
        values.append(df_genes[i])

    return go.Figure(
        data=[
            go.Table(
                header={
                    "values": c_name,
                    "fill_color": "paleturquoise",
                    "align": "left",
                },
                cells={
                    "values": values,
                    "fill_color": "lavender",
                    "align": "left",
                },
            ),
        ],
    )


# MULTI IMAGE BOX/PIE
def func_multi_plot(list_clusters: list, col_name: str):
    fig = None
    df_clinical = DF_CLINICAL_DATA[
        DF_CLINICAL_DATA["cluster"].isin(list_clusters)
    ]
    df_clinical = df_clinical.dropna(subset=[col_name])
    df_clinical["cluster_plot"] = df_clinical["cluster"].apply(
        lambda x: f"cluster_{x}",
    )

    if col_name in NUMERIC_COLUMNS_CLINICAL:
        fig = tap.plot_stats(df_clinical, "cluster_plot", col_name)
    else:
        fig = make_subplots(
            1,
            len(list_clusters),
            subplot_titles=[f"cluster {e}" for e in list_clusters],
            specs=[[{"type": "domain"} for e in list_clusters]],
        )
        for i, index in enumerate(list_clusters):
            _cluster_values = df_clinical[df_clinical["cluster"] == index][
                col_name
            ]
            _temp_dict = dict(_cluster_values.value_counts())
            fig.add_trace(
                go.Pie(
                    labels=list(_temp_dict.keys()),
                    values=list(_temp_dict.values()),
                    scalegroup="one",
                ),
                1,
                i + 1,
            )
    return fig


# UPDATE MULTI FIG1
@callback(
    Output("fig_multi_fig", "figure"),
    [Input("dd-cluster-multi", "value"), Input("dd-box-1", "value")],
)
def update_multi_fig(list_clusters: list, col_name: str):
    global CLUSTER_SELECTED_MULTI
    global BOX_FIG_SELECTED_1
    min_cluster = 2
    if len(list_clusters) < min_cluster or col_name is None:
        return no_update
    CLUSTER_SELECTED_MULTI = list_clusters
    BOX_FIG_SELECTED_1 = col_name
    return func_multi_plot(list_clusters, col_name)


# START
if __name__ == "__main__":
    APP.run(debug=False, host="0.0.0.0", port=8593)
