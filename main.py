"""
This script sets up a Dash web application for visualizing and analyzing data.
The application includes multiple pages for different types of analysis,
including study creation, study description, pathway analysis,
clinical data visualization, survival analysis, and cluster comparison.
Modules:
    - base64
    - io
    - json
    - os
    - pickle
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
    - redirect_pages(pathname): Redirects to the appropriate page based on the URL pathname.
    - update_dropdown_liststudy(n_clicks): Updates the dropdown list of studies.
    - select_study(value): Selects a study and loads its data.
    - update_mutational_file(filename): Updates the mutational file upload component.
    - update_clinical_patient_file(filename): Updates the clinical patient file upload component.
    - update_clinical_sample_file(filename): Updates the clinical sample file upload component.
    - update_val_col_patient_name(data, sep, skip): Updates the dropdown options.
    - update_val_col_sample_name(data, sep, skip): Updates the dropdown options.
    - update_list_columns_mutation(data, sep, skip): Updates the dropdown options.
    - create_study(n_clicks, ...): Creates a new study based on the provided input data.
    - dropdpwn_cluster(n_clicks): Updates the cluster dropdown options.
    - displaySelectedNodeData(data_dict): Displays data for the selected node in the graph.
    - update_graph(layout): Updates the graph layout.
    - update_cluster(cluster): Updates the cluster elements and related figures.
    - update_go(cluster, pvalue, adjusted_pvalue, process_type): Updates the GO figure.
    - update_kegg(cluster, pvalue, adjusted_pvalue): Updates the KEGG figure.
    - update_reactome(cluster, pvalue, adjusted_pvalue): Updates the REACTOME figure.
    - update_wiki(cluster, pvalue, adjusted_pvalue): Updates the WIKI figure.
    - dropdown_box_fig_1(n_clicks): Updates the dropdown options.
    - dropdown_box_fig_2(n_clicks): Updates the dropdown options.
    - func_single_plot(cluster, column_name): Generates a single plot for the cluster.
    - update_box_1(cluster, column_name): Updates the first box plot.
    - update_box_2(cluster, column_name): Updates the second box plot.
    - update_table_clinical_data(cluster): Updates the clinical data table.
    - update_overall_survival(cluster): Updates the overall survival plot and statistics.
    - update_survival_comparison(list_clusters): Updates the survival comparison plot.
    - update_venn(list_clusters): Updates the Venn diagram for gene comparison.
    - update_genes_common(list_clusters): Updates the table of common genes.
    - func_multi_plot(list_clusters, column_name): Generates a multi-cluster plot for the column.
    - update_multi_fig1(list_clusters, column_name): Updates the first multi-cluster plot.
Usage:
    Run this script to start the Dash web application.
    The application will be available at http://127.0.0.1:8593.
"""

import base64
import io
import json
import os
import pickle

import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import tap
from dash import Dash, Input, Output, State, callback
from dash import dash_table, dcc, html, no_update
from lifelines import KaplanMeierFitter
from lifelines.statistics import pairwise_logrank_test
from plotly.subplots import make_subplots

from lib import venn

# Config lib
matplotlib.use("agg")
pd.options.mode.copy_on_write = True
cyto.load_extra_layouts()
# Global Vars
APP_NAME = "Tortoise"
APP_LOGO = "/assets/tortoise.png"
LIST_STUDIES = []
NAME_STUDY = None
STUDY_CONFIG = None
OUT_R_PATH = ""
DF_CLINICAL_DATA = None
NUMERIC_COLUMNS_CLINICAL = []
ALL_COLUMNS_CLINICAL = []
GRAPH = None
GRAPH_ELEMENTS = None
CLUSTERS_INDEX = []
GENES = []
CLUSTER_SELECTED = None
CLUSTER_SELECTED_MULTI = []
BOX_FIG_SELECTED_1 = None
BOX_FIG_SELECTED_2 = None
ALL_GENES = []


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
            }
        )
    # convert and add edges
    for e in graph_filtered.es():
        g_ele.append({"data": {"source": e.source, "target": e.target}})
    return g_ele


# APP + SIDEBAR
APP = Dash(
    title=APP_NAME, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME]
)
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
                        html.Span("Study Description", className="navbar_span"),
                    ],
                    href="/study_description",
                    active="exact",
                    className="navbar_entity",
                ),
                # PATWAY
                dbc.NavLink(
                    [
                        html.I(className="fas fa-diagram-project"),
                        html.Span("Pathways Analysis", className="navbar_span"),
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
                        html.Span("Survival Analysis", className="navbar_span"),
                    ],
                    href="/survival_analysis",
                    active="exact",
                    className="navbar_entity",
                ),
                # CLUSTER COMP
                dbc.NavLink(
                    [
                        html.I(className="fas fa-code-compare"),
                        html.Span("Cluster Comparision", className="navbar_span"),
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
    dbc.Row([html.H2(f"Welcome to {APP_NAME}", style={"text-align": "center"})]),
    html.Hr(),
    dbc.Row(
        [
            dbc.Col([html.H5("Selected Study")], lg=2),
            dbc.Col([dcc.Dropdown(id="dd-study")], lg=4),
        ],
        justify="center",
    ),
    html.Hr(),
    dbc.Row(
        [
            html.Span(
                id="dd-study-info",
                className="dd_study_info",
                style={"text-align": "center"},
            )
        ]
    ),
]


# PAGING
@callback(Output("page-content", "children"), Input("url", "pathname"))
def redirect_pages(pathname):
    if pathname == "/":
        return PAGE_HOME
    elif pathname == "/create_study":
        return PAGE_CREATE_STUDY
    elif pathname == "/study_description":
        return PAGE_STUDY_DESCRIPTION
    elif pathname == "/pathway_analysis":
        return PAGE_PATHWAY_ANALYSIS
    elif pathname == "/cluster_comparision":
        return PAGE_CLUSTER_COMPARISION
    elif pathname == "/clinical_data":
        return PAGE_CLINICAL_DATA
    elif pathname == "/survival_analysis":
        return PAGE_SURVIVAL_ANALYSIS
    else:
        # If the user tries to reach a different page, return a 404 message
        return html.Div(
            [
                html.H1("404: Not found", className="text-danger"),
                html.Hr(),
                html.P(f"The pathname {pathname} was not recognised..."),
            ],
            className="p-3 bg-light rounded-3",
        )


@callback(
    Output("dd-study", "options"),
    Output("dd-study", "value"),
    Input("dd-study", "n_clicks"),
)
def update_dropdown_liststudy(n_clicks):
    global LIST_STUDIES
    if not os.path.exists("study"):
        LIST_STUDIES = []
    else:
        LIST_STUDIES = os.listdir("study")
    return LIST_STUDIES, NAME_STUDY


# SELECT STUDY
@callback(
    Output("dd-study-info", "children"),
    Input("dd-study", "value"),
    prevent_initial_call=True,
)
def select_study(value):
    global NAME_STUDY
    global OUT_R_PATH
    global GRAPH
    global CLUSTERS_INDEX
    global DF_CLINICAL_DATA
    global NUMERIC_COLUMNS_CLINICAL
    global ALL_COLUMNS_CLINICAL
    global CLUSTER_SELECTED
    global CLUSTER_SELECTED_MULTI
    global BOX_FIG_SELECTED_1
    global BOX_FIG_SELECTED_2
    global STUDY_CONFIG
    global ALL_GENES

    if value is None:
        return ""
    # prevent reload
    if value == NAME_STUDY:
        return ""

    NAME_STUDY = value
    OUT_R_PATH = os.path.join("study", NAME_STUDY, "output")
    with open(os.path.join(OUT_R_PATH, "graph.pickle"), "rb") as f:
        GRAPH = pickle.load(f)

    CLUSTERS_INDEX = [int(c) for c in set(GRAPH.vs["cluster"])]
    ALL_GENES = set(GRAPH.vs["gene"])

    # Find all numerical column
    DF_CLINICAL_DATA = pd.read_csv(
        os.path.join(OUT_R_PATH, "cluster_clinical_data.csv"),
        sep="\t",
        engine="python",
    )
    DF_CLINICAL_DATA_ALL = DF_CLINICAL_DATA.copy()
    DF_CLINICAL_DATA_ALL["cluster"] = "ALL"
    DF_CLINICAL_DATA = pd.concat([DF_CLINICAL_DATA, DF_CLINICAL_DATA_ALL])
    NUMERIC_COLUMNS_CLINICAL = list(
        DF_CLINICAL_DATA.select_dtypes(include=np.number).columns
    )
    ALL_COLUMNS_CLINICAL = list(DF_CLINICAL_DATA.columns)
    ALL_COLUMNS_CLINICAL.remove("cluster")
    # Load Config
    with open(os.path.join("study", NAME_STUDY, "config.json"), "r") as config_file:
        STUDY_CONFIG = json.load(config_file)
    # Reset selection
    CLUSTER_SELECTED = CLUSTERS_INDEX[0]
    CLUSTER_SELECTED_MULTI = ["ALL", CLUSTERS_INDEX[0]]
    BOX_FIG_SELECTED_1 = ALL_COLUMNS_CLINICAL[0]
    BOX_FIG_SELECTED_2 = ALL_COLUMNS_CLINICAL[-1]
    return f'You have selected "{value}" study'


# CREATE STUDY
PAGE_CREATE_STUDY = [
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
                        id="input_namestudy", type="text", placeholder="Name study"
                    ),
                ],
                width=2,
            ),
        ]
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
                    dcc.Upload(
                        id="mutational-file",
                        children=html.Span("Upload Mutational File"),
                        # Allow multiple files to be uploaded
                        multiple=False,
                    )
                ],
                width=2,
                className="upload-button",
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
        ]
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
                    )
                ],
                width=2,
            ),
        ]
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
                    )
                ],
                width=2,
            ),
        ]
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
                [dcc.Dropdown(options=[], id="dd-identifier-columns", multi=True)],
                width=2,
            ),
        ]
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
                    )
                ],
                width=2,
            ),
        ]
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
                    ),
                ],
                width=4,
            ),
        ]
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
                    dcc.Upload(
                        id="clinical-patient-file",
                        children=html.Span("Upload Clinical Patient File"),
                        # Allow multiple files to be uploaded
                        multiple=False,
                    )
                ],
                width=2,
                className="upload-button",
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
        ]
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
                    )
                ],
                width=2,
            ),
        ]
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
                    dcc.Upload(
                        id="clinical-sample-file",
                        children=html.Span("Upload Clinical Sample File"),
                        # Allow multiple files to be uploaded
                        multiple=False,
                    )
                ],
                width=2,
                className="upload-button",
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
        ]
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
                    )
                ],
                width=2,
            ),
        ]
    ),
    # CREATE STUDY BUTTON
    dbc.Row([html.Button("Create Study", id="create-study-button", n_clicks=0)]),
    # LOADING
    dbc.Row(
        [
            dcc.Loading(
                children=[
                    html.Span(id="create-study-log", className="create_study_log")
                ],
            ),
        ]
    ),
]


@callback(
    Output("mutational-file", "children"),
    Input("mutational-file", "filename"),
    prevent_initial_call=True,
)
def update_mutational_file(filename):
    return filename


@callback(
    Output("clinical-patient-file", "children"),
    Input("clinical-patient-file", "filename"),
    prevent_initial_call=True,
)
def update_clinical_patient_file(filename):
    return filename


@callback(
    Output("clinical-sample-file", "children"),
    Input("clinical-sample-file", "filename"),
    prevent_initial_call=True,
)
def update_clinical_sample_file(filename):
    return filename


# dropdown column patient name
@callback(
    Output("dd-column-patient-name", "options"),
    [
        Input("clinical-patient-file", "contents"),
        Input("clinical-patient-separator", "value"),
        Input("clinical-patient-skiprow", "value"),
    ],
    prevent_initial_call=True,
)
def update_val_col_patient_name(data, sep, skip):
    if data is None or sep is None or skip is None:
        return []
    content_type, content_string = data.split(",")
    decoded = base64.b64decode(content_string)
    df = pd.read_csv(
        io.StringIO(decoded.decode("utf-8")), sep=sep, skiprows=skip, engine="python"
    )
    return df.columns


# dropdown column sample name
@callback(
    Output("dd-column-sample-name", "options"),
    [
        Input("clinical-sample-file", "contents"),
        Input("clinical-sample-separator", "value"),
        Input("clinical-sample-skiprow", "value"),
    ],
    prevent_initial_call=True,
)
def update_val_col_sample_name(data, sep, skip):
    if data is None or sep is None or skip is None:
        return []
    content_type, content_string = data.split(",")
    decoded = base64.b64decode(content_string)
    df = pd.read_csv(
        io.StringIO(decoded.decode("utf-8")), sep=sep, skiprows=skip, engine="python"
    )
    return df.columns


# dropdown mutation
@callback(
    Output("dd-column-gene", "options"),
    Output("dd-identifier-columns", "options"),
    Output("dd-vaf", "options"),
    Output("dd-column-sample-name-mutation", "options"),
    [
        Input("mutational-file", "contents"),
        Input("mutational-separator", "value"),
        Input("mutational-skiprow", "value"),
    ],
    prevent_initial_call=True,
)
def update_list_columns_mutation(data, sep, skip):
    if data is None or sep is None or skip is None:
        return [], [], [], []
    content_type, content_string = data.split(",")
    decoded = base64.b64decode(content_string)
    df = pd.read_csv(
        io.StringIO(decoded.decode("utf-8")), sep=sep, skiprows=skip, engine="python"
    )
    return df.columns, df.columns, df.columns, df.columns


@callback(
    Output("create-study-log", "children"),
    Input("create-study-button", "n_clicks"),
    State("input_namestudy", "value"),
    State("mutational-file", "contents"),
    State("mutational-separator", "value"),
    State("mutational-skiprow", "value"),
    State("clinical-patient-file", "contents"),
    State("clinical-patient-separator", "value"),
    State("clinical-patient-skiprow", "value"),
    State("clinical-sample-file", "contents"),
    State("clinical-sample-separator", "value"),
    State("clinical-sample-skiprow", "value"),
    State("dd-column-patient-name", "value"),
    State("dd-column-sample-name", "value"),
    State("dd-column-sample-name-mutation", "value"),
    State("dd-column-gene", "value"),
    State("dd-identifier-columns", "value"),
    State("dd-vaf", "value"),
    State("input-vaf-score", "value"),
    prevent_initial_call=True,
)
def create_study(
    n_clicks,
    input_namestudy,
    mutational_data,
    mutational_separator,
    mutational_skiprow,
    clinical_patient_data,
    clinical_patient_separator,
    clinical_patient_skiprow,
    clinical_sample_data,
    clinical_sample_separator,
    clinical_sample_skiprow,
    c_patient_name,
    c_sample_name,
    c_sample_mutation,
    c_gene,
    c_identifier,
    c_vaf,
    vaf_score,
):
    # CHECK INPUT DATA
    if input_namestudy is None:
        return html.Div("Please insert name for study")
    if os.path.exists(os.path.join("study", input_namestudy, "input")):
        return html.Div("Name study already exist")
    if mutational_data is None:
        return html.Div("Please insert file for mutational data")
    if mutational_separator is None:
        return html.Div("Please select separator for mutational data")
    if mutational_skiprow is None:
        return html.Div("Please select skiprow for mutational data")
    if c_gene is None:
        return html.Div("Please select column for gene")
    if c_sample_mutation is None:
        return html.Div("Please select column for sample name on mutation file")
    if c_identifier == []:
        return html.Div("Please select at least one column for mutation identifier")
    if clinical_patient_data is not None:
        if clinical_patient_separator is None:
            return html.Div("Please select separator for clinical patient data")
        if clinical_patient_skiprow is None:
            return html.Div("Please select skiprow for clinical patient data")
        if c_patient_name is None:
            return html.Div("Please select column for patient name")
    if clinical_sample_data is not None:
        if clinical_sample_separator is None:
            return html.Div("Please select separator for clinical sample data")
        if clinical_sample_skiprow is None:
            return html.Div("Please select skiprow for clinical sample data")
        if c_sample_name is None:
            return html.Div(
                "Please select column for sample name on clinical sample file"
            )
    # CREATE STUDY FOLDER
    os.makedirs(os.path.join("study", input_namestudy, "input"))
    # MUTATIONAL DATA
    content_type, content_string = mutational_data.split(",")
    decoded = base64.b64decode(content_string)
    df = pd.read_csv(
        io.StringIO(decoded.decode("utf-8")),
        sep=mutational_separator,
        skiprows=mutational_skiprow,
        engine="python",
    )
    df.to_csv(
        os.path.join("study", input_namestudy, "input", "mutational_data.txt"),
        sep="\t",
        index=False,
    )
    # CLINICAL PATIENT
    if clinical_patient_data is not None:
        content_type, content_string = clinical_patient_data.split(",")
        decoded = base64.b64decode(content_string)
        df = pd.read_csv(
            io.StringIO(decoded.decode("utf-8")),
            sep=clinical_patient_separator,
            skiprows=clinical_patient_skiprow,
            engine="python",
        )
        df.to_csv(
            os.path.join("study", input_namestudy, "input", "clinical_patient.txt"),
            sep="\t",
            index=False,
        )
    # CLINICAL SAMPLE
    if clinical_sample_data is not None:
        content_type, content_string = clinical_sample_data.split(",")
        decoded = base64.b64decode(content_string)
        df = pd.read_csv(
            io.StringIO(decoded.decode("utf-8")),
            sep=clinical_sample_separator,
            skiprows=clinical_sample_skiprow,
            engine="python",
        )
        df.to_csv(
            os.path.join("study", input_namestudy, "input", "clinical_sample.txt"),
            sep="\t",
            index=False,
        )
    # GENERATE JSON CONFIG
    DATA = {}
    DATA["Paths"] = {}
    DATA["Clinical_data"] = {}
    DATA["Mutation"] = {}

    DATA["Paths"]["name_study"] = input_namestudy
    DATA["Paths"]["data_mutational"] = os.path.join(
        "study", input_namestudy, "input", "mutational_data.txt"
    )
    DATA["Paths"]["data_mutational_sep"] = "\t"
    DATA["Paths"]["data_mutational_skip"] = 0
    DATA["Paths"]["data_clinical_patient"] = ""
    DATA["Paths"]["data_clinical_sample_sep"] = "\t"
    DATA["Paths"]["data_clinical_sample_skip"] = 0
    DATA["Paths"]["data_clinical_sample"] = ""
    DATA["Paths"]["data_clinical_patient_sep"] = "\t"
    DATA["Paths"]["data_clinical_patient_skip"] = 0
    DATA["Mutation"]["column_gene"] = c_gene
    DATA["Mutation"]["column_sample_name"] = c_sample_mutation
    # REMOVE EXTRA SEPARATOR BEFOR JOIN
    c_identifier_clean = [col.replace(";", "") for col in c_identifier]
    DATA["Mutation"]["identifier_columns"] = ";".join(c_identifier_clean)
    DATA["Clinical_data"]["column_patient_name"] = ""
    DATA["Clinical_data"]["column_sample_name"] = ""
    DATA["Mutation"]["vaf"] = False if vaf_score is None else True
    DATA["Mutation"]["vaf_score"] = vaf_score
    DATA["Mutation"]["vaf_column"] = "" if c_vaf is None else c_vaf

    if clinical_patient_data is not None:
        DATA["Paths"]["data_clinical_patient"] = os.path.join(
            "study", input_namestudy, "input", "clinical_patient.txt"
        )
        DATA["Clinical_data"]["column_patient_name"] = c_patient_name
    if clinical_sample_data is not None:
        DATA["Paths"]["data_clinical_sample"] = os.path.join(
            "study", input_namestudy, "input", "clinical_sample.txt"
        )
        DATA["Clinical_data"]["column_sample_name"] = c_sample_name
    path_config = os.path.join("study", input_namestudy, "config.json")
    with open(path_config, "w") as f:
        json.dump(DATA, f, indent=4)

    print(os.popen(f"python tortoise.py -c {path_config}").read())
    return html.Div(f"{input_namestudy} -- Study created!")


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
        ]
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
                                "style": {"content": "data(name)", "font-size": "5px"},
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
                            "concentric",
                            "cose",
                            "grid",
                            "circle",
                            "breadthfirst",
                            "klay",
                        ],  # cola
                        value="concentric",
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
                    )
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
                    )
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
def dropdpwn_cluster(n_clicks):
    global CLUSTERS_INDEX
    global CLUSTER_SELECTED
    return CLUSTERS_INDEX, CLUSTER_SELECTED


# SELECT SINGLE NODE GRAPH
@callback(
    Output("modal-lg", "children"),
    Output("modal-lg", "is_open"),
    Input("cytoscape-graph", "tapNodeData"),
    prevent_initial_call=True,
)
def displaySelectedNodeData(data_dict):
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
        os.path.join(OUT_R_PATH, "Gene_Count", f"genes_cluster_{cluster}.csv"),
        sep="\t",
        engine="python",
    )
    fig_pie = px.pie(
        df_gene, values="COUNT", names="GENE", title="Number Mutation for Gene"
    )
    fig_pie.update_traces(textposition="inside", textinfo="label")
    # VARIANT NUMBERS
    df_variant = pd.read_csv(
        os.path.join(
            OUT_R_PATH, "Variants_Degree", f"variants_degree_cluster{cluster}.csv"
        ),
        sep="\t",
        engine="python",
    )
    n_variants = len(df_variant)
    # FIGURE DEGREE
    df_variant = df_variant.sort_values(by=["Degree"], ascending=False)[:15]
    fig_degree = px.bar(df_variant, x="Variants", y="Degree", title="Mutation Degree")
    # PATIENTS NUMBER
    n_patients = len(
        [1 for e in cluster_elements if e["data"].get("vertex_type", "") == "PATIENT"]
    )
    # GENE NUMBERS
    n_genes = len(df_gene)
    # VARIANT CENTROID
    if n_variants > 1 and df_variant.iloc[0]["Degree"] == df_variant.iloc[1]["Degree"]:
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
                    html.Span("Use adjusted PValue", className="span_selector"),
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
        ]
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
                            {"label": "Biological Function", "value": "biological"},
                            {"label": "Molecular Function", "value": "molecular"},
                            {"label": "Cellular Component", "value": "cellular"},
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
                    )
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
                    )
                ],
                lg=6,
            ),
        ]
    ),
]


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
def update_go(cluster, pvalue, adjusted_pvalue, process_type):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df = pd.read_csv(
        os.path.join(
            OUT_R_PATH,
            "Arricchimento_all_genes",
            "GO",
            f"{process_type}_{cluster}.csv",
        ),
        engine="python",
    )
    if adjusted_pvalue == "True":
        df = df[df["Adjusted.P.value"] < pvalue]
        df = df.sort_values(by=["Adjusted.P.value"], ascending=False)[-25:]
        df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
        fig = px.bar(
            df,
            x="Count_gene",
            y="Term",
            hover_data=["Overlap"],
            color="Adjusted.P.value",
            title="GO",
            labels={"Adjusted.P.value": "Adjusted Pvalue"},
        )
        fig.update_layout(
            xaxis_title="Genes_Count",
            yaxis_title="Terms",
            legend_title="Adjusted Pvalue",
            coloraxis_colorbar={"exponentformat": "e"},
        )
        return fig

    df = df[df["P.value"] < pvalue]
    df = df.sort_values(by=["P.value"], ascending=False)[-25:]
    df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
    fig = px.bar(
        df,
        x="Count_gene",
        y="Term",
        hover_data=["Overlap"],
        color="P.value",
        title="GO",
        labels={"P.value": "Pvalue"},
    )
    fig.update_layout(
        xaxis_title="Genes_Count",
        yaxis_title="Terms",
        legend_title="Pvalue",
        coloraxis_colorbar={"exponentformat": "e"},
    )
    return fig


# UPDATE KEGG FIGURE
@callback(
    Output("fig_kegg", "figure"),
    [
        Input("dd-cluster", "value"),
        Input("dd-pvalue", "value"),
        Input("dd-adjusted-pvalue", "value"),
    ],
)
def update_kegg(cluster, pvalue, adjusted_pvalue):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df = pd.read_csv(
        os.path.join(
            OUT_R_PATH, "Arricchimento_all_genes", "KEGG", f"kegg_{cluster}.csv"
        ),
        engine="python",
    )
    if adjusted_pvalue == "True":
        df = df[df["Adjusted.P.value"] < pvalue]
        df = df.sort_values(by=["Adjusted.P.value"], ascending=False)[-25:]
        df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
        fig = px.bar(
            df,
            x="Count_gene",
            y="Term",
            hover_data=["Overlap"],
            color="Adjusted.P.value",
            title="KEGG",
            color_continuous_scale=px.colors.sequential.Viridis,
            labels={"P.value": "Adjusted Pvalue"},
        )
        fig.update_layout(
            xaxis_title="Genes_Count",  # Nome dell'asse delle x
            yaxis_title="Terms",
            legend_title="Adjusted Pvalue",
            coloraxis_colorbar={"exponentformat": "e"},
        )
        return fig
    else:
        df = df[df["P.value"] < pvalue]
        df = df.sort_values(by=["P.value"], ascending=False)[-25:]
        df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
        fig = px.bar(
            df,
            x="Count_gene",
            y="Term",
            hover_data=["Overlap"],
            color="P.value",
            title="KEGG",
            color_continuous_scale=px.colors.sequential.Viridis,
            labels={"P.value": "Pvalue"},
        )
        fig.update_layout(
            xaxis_title="Genes_Count",  # Nome dell'asse delle x
            yaxis_title="Terms",
            legend_title="Pvalue",
            coloraxis_colorbar={"exponentformat": "e"},
        )
        return fig


# UPDATE REACTOME FIGURE
@callback(
    Output("fig_reactome", "figure"),
    [
        Input("dd-cluster", "value"),
        Input("dd-pvalue", "value"),
        Input("dd-adjusted-pvalue", "value"),
    ],
)
def update_reactome(cluster, pvalue, adjusted_pvalue):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df = pd.read_csv(
        os.path.join(
            OUT_R_PATH,
            "Arricchimento_all_genes",
            "REACTOME",
            f"reactome_{cluster}.csv",
        ),
        engine="python",
    )
    if adjusted_pvalue == "True":
        df = df[df["Adjusted.P.value"] < pvalue]
        df = df.sort_values(by=["Adjusted.P.value"], ascending=False)[-25:]
        df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
        fig = px.bar(
            df,
            x="Count_gene",
            y="Term",
            hover_data=["Overlap"],
            color="Adjusted.P.value",
            title="REACTOME",
            color_continuous_scale=px.colors.sequential.Viridis,
            labels={"P.value": "Adjusted Pvalue"},
        )
        fig.update_layout(
            xaxis_title="Genes_Count",  # Nome dell'asse delle x
            yaxis_title="Terms",
            legend_title="Adjusted Pvalue",
            coloraxis_colorbar={"exponentformat": "e"},
        )
        return fig
    else:
        df = df[df["P.value"] < pvalue]
        df = df.sort_values(by=["P.value"], ascending=False)[-25:]
        df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
        fig = px.bar(
            df,
            x="Count_gene",
            y="Term",
            hover_data=["Overlap"],
            color="P.value",
            title="REACTOME",
            color_continuous_scale=px.colors.sequential.Viridis,
            labels={"P.value": "Pvalue"},
        )
        fig.update_layout(
            xaxis_title="Genes_Count",  # Nome dell'asse delle x
            yaxis_title="Terms",
            legend_title="Pvalue",
            coloraxis_colorbar={"exponentformat": "e"},
        )
        return fig


# UPDATE WIKI FIGURE
@callback(
    Output("fig_wiki", "figure"),
    [
        Input("dd-cluster", "value"),
        Input("dd-pvalue", "value"),
        Input("dd-adjusted-pvalue", "value"),
    ],
)
def update_wiki(cluster, pvalue, adjusted_pvalue):
    if cluster is None:
        return no_update
    global CLUSTER_SELECTED
    CLUSTER_SELECTED = cluster
    df = pd.read_csv(
        os.path.join(
            OUT_R_PATH, "Arricchimento_all_genes", "WIKI", f"wiki_{cluster}.csv"
        ),
        engine="python",
    )
    if adjusted_pvalue == "True":
        df = df[df["Adjusted.P.value"] < pvalue]
        df = df.sort_values(by=["Adjusted.P.value"], ascending=False)[-25:]
        df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
        fig = px.bar(
            df,
            x="Count_gene",
            y="Term",
            hover_data=["Overlap"],
            color="Adjusted.P.value",
            title="WIKI",
            color_continuous_scale=px.colors.sequential.Viridis,
            labels={"Adjusted.P.value": "Adjusted Pvalue", "Overlap": "Overlap_Genes"},
        )
        fig.update_layout(
            xaxis_title="Genes_Count",  # Nome dell'asse delle x
            yaxis_title="Terms",
            legend_title="Adjusted Pvalue",
            coloraxis_colorbar={"exponentformat": "e"},
        )
        return fig
    else:
        df = df[df["P.value"] < pvalue]
        df = df.sort_values(by=["P.value"], ascending=False)[-25:]
        df["Count_gene"] = df.apply(lambda e: len(e["Genes"].split(";")), axis=1)
        fig = px.bar(
            df,
            x="Count_gene",
            y="Term",
            hover_data=["Overlap"],
            color="P.value",
            title="WIKI",
            color_continuous_scale=px.colors.sequential.Viridis,
            labels={"P.value": "Pvalue", "Overlap": "Overlap_Genes"},
        )
        fig.update_layout(
            xaxis_title="Genes_Count",  # Nome dell'asse delle x
            yaxis_title="Terms",
            legend_title="Pvalue",
            coloraxis_colorbar={"exponentformat": "e"},
        )
        return fig


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
        ]
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
        ]
    ),
    dbc.Row([dash_table.DataTable(id="table_clinical_data")]),
]


@callback(
    Output("dd-box-1", "options"),
    Output("dd-box-1", "value"),
    Input("dd-box-1", "n_clicks"),
)
def dropdown_box_fig_1(n_clicks):
    global ALL_COLUMNS_CLINICAL
    global BOX_FIG_SELECTED_1
    return ALL_COLUMNS_CLINICAL, BOX_FIG_SELECTED_1


@callback(
    Output("dd-box-2", "options"),
    Output("dd-box-2", "value"),
    Input("dd-box-2", "n_clicks"),
)
def dropdown_box_fig_2(n_clicks):
    global ALL_COLUMNS_CLINICAL
    global BOX_FIG_SELECTED_2
    return ALL_COLUMNS_CLINICAL, BOX_FIG_SELECTED_2


# SINGLE IMAGE BOX/PIE
def func_single_plot(cluster, column_name):
    cluster_values = DF_CLINICAL_DATA[DF_CLINICAL_DATA["cluster"] == cluster][
        column_name
    ]

    if column_name in NUMERIC_COLUMNS_CLINICAL:
        fig = px.box(cluster_values, y=column_name)
        return fig
    else:
        _temp_dict = dict(cluster_values.value_counts())
        _temp_df = pd.DataFrame(
            {column_name: _temp_dict.keys(), "count": _temp_dict.values()}
        )
        fig = px.pie(_temp_df, values="count", names=column_name)
        return fig


# UPDATE BOX_PLOT_1
@callback(
    Output("fig_box_plot_1", "figure"),
    [Input("dd-cluster", "value"), Input("dd-box-1", "value")],
)
def update_box_1(cluster, column_name):
    if cluster is None or column_name is None:
        return no_update
    global BOX_FIG_SELECTED_1
    global CLUSTER_SELECTED
    BOX_FIG_SELECTED_1 = column_name
    CLUSTER_SELECTED = cluster
    return func_single_plot(cluster, column_name)


# UPDATE BOX_PLOT_2
@callback(
    Output("fig_box_plot_2", "figure"),
    [Input("dd-cluster", "value"), Input("dd-box-2", "value")],
)
def update_box_2(cluster, column_name):
    if cluster is None or column_name is None:
        return no_update
    global BOX_FIG_SELECTED_2
    global CLUSTER_SELECTED
    BOX_FIG_SELECTED_2 = column_name
    CLUSTER_SELECTED = cluster
    return func_single_plot(cluster, column_name)


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
        ]
    ),
    dbc.Row(
        [
            dbc.Col(
                [
                    dcc.Graph(
                        id="survival_figure",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    )
                ],
                lg=8,
            ),
            dbc.Col(
                [
                    dcc.Graph(
                        id="survival_figure_stat",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    )
                ],
                lg=4,
            ),
        ]
    ),
    dbc.Row(
        [
            # CLUSTER SELECTOR MULTI
            dbc.Col(
                [
                    html.Span("Multi Cluster Selection", className="span_selector"),
                    dcc.Dropdown(
                        id="dd-cluster-multi",
                        multi=True,
                    ),
                    html.Br(),
                ],
                lg=4,
            ),
            html.Hr(),
        ]
    ),
    dbc.Row(
        [
            dbc.Col(
                [
                    dcc.Graph(
                        id="survival_figure_comparison",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    )
                ],
                lg=8,
            ),
            dbc.Col(
                [
                    dcc.Graph(
                        id="table_test_survival",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    )
                ],
                lg=4,
            ),
        ]
    ),
]


@callback(
    Output("dd-cluster-multi", "options"),
    Output("dd-cluster-multi", "value"),
    Input("dd-cluster-multi", "n_clicks"),
)
def dropdpwn_multi_cluster(n_clicks):
    global CLUSTERS_INDEX
    global CLUSTER_SELECTED_MULTI
    return ["ALL"] + CLUSTERS_INDEX, CLUSTER_SELECTED_MULTI


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
    # FIXME hardcoded name column
    col_status = "VITAL_STATUS"
    column_name_month = "OS_INT"
    data = DF_CLINICAL_DATA[DF_CLINICAL_DATA["cluster"] == cluster]
    data.dropna(subset=[column_name_month, col_status], inplace=True)
    if len(data) <= 0:
        return no_update, no_update
    data[col_status] = data[col_status].replace({"Yes": 1, "No": 0})
    kmf = KaplanMeierFitter()
    kmf.fit(data[column_name_month].values, event_observed=data[col_status].values)
    # Crea il grafico della curva di sopravvivenza con Plotly
    fig = go.Figure()
    # Aggiungi la curva di sopravvivenza
    fig.add_trace(
        go.Scatter(
            x=kmf.confidence_interval_.index,
            y=kmf.confidence_interval_["KM_estimate_upper_0.95"],
            mode="lines",
            line={"shape": "hv", "width": 0},
            showlegend=False,
        )
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
        )
    )
    fig.update_layout(
        title=f"Kaplan-Meier Curve Cluster: {cluster}",
        xaxis_title="Duration",
        yaxis_title="Survival probability",
        # margin=dict(r=0, t=10, l=0),
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
        )
    )
    # STAT PIE
    fig_stat = func_single_plot(cluster, col_status)
    fig_stat.update_layout(
        {"title": f"Vital Status Cluster: {cluster}", "legend_title": "Alive"}
    )
    return fig, fig_stat


@callback(
    Output("survival_figure_comparison", "figure"),
    Output("table_test_survival", "figure"),
    Input("dd-cluster-multi", "value"),
)
def update_survival_comparison(list_clusters):
    global CLUSTER_SELECTED_MULTI
    if len(list_clusters) < 2:
        return no_update, no_update
    CLUSTER_SELECTED_MULTI = list_clusters
    fig = go.Figure()
    fig_stats = None
    kmf = KaplanMeierFitter()
    # FIXME hardcoded name column
    col_status = "VITAL_STATUS"
    column_name_month = "OS_INT"

    DF = DF_CLINICAL_DATA[DF_CLINICAL_DATA["cluster"].isin(list_clusters)]
    DF.dropna(subset=[column_name_month, col_status], inplace=True)
    DF[col_status] = DF[col_status].replace({"Yes": 1, "No": 0})

    for cluster in list(DF["cluster"].unique()):
        cluster_data = DF[DF["cluster"] == cluster]
        if len(cluster_data) == 0:
            list_clusters.remove(cluster)
        else:
            # Fit e plot della curva di sopravvivenza
            kmf.fit(
                cluster_data[column_name_month],
                event_observed=cluster_data[col_status],
                label=f"Cluster {cluster}",
            )
            # kmf.plot_survival_function(ci_show=True)
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
                )
            )

    if len(list_clusters) <= 1:
        return no_update, no_update

    DF["cluster_str"] = DF["cluster"].astype(str)
    # Test log-rank per confrontare le curve tra tutti i cluster
    results = pairwise_logrank_test(
        DF[column_name_month], DF["cluster_str"], DF[col_status]
    )
    df_res = results.summary
    df_res = df_res.rename(columns={"-log2(p)": "log2_p"})
    df_res = df_res.reset_index()
    df_res["comparison"] = (
        df_res["level_0"].astype(str) + " " + df_res["level_1"].astype(str)
    )
    # Rimuovi le colonne level_0 e level_1, se non necessarie
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
            )
        ]
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
        ]
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
                    )
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
                        id="fig_multi_fig1",
                        className="add-border",
                        style={"width": "100%", "height": "38vh"},
                    ),
                ],
                lg=6,
            ),
        ]
    ),
]


# UPDATE COMPARISION
@callback(Output("plot-venn", "src"), Input("dd-cluster-multi", "value"))
def update_venn(list_clusters):
    global CLUSTER_SELECTED_MULTI
    if len(list_clusters) < 2:
        return no_update
    CLUSTER_SELECTED_MULTI = list_clusters

    cluster_gene_list = []
    for index in list_clusters:
        gene_values = []
        if index == "ALL":
            gene_values = pd.read_csv(
                os.path.join(OUT_R_PATH, "distribution_gene_cluster.csv"),
                sep="\t",
                engine="python",
            )["Gene"].unique()
        else:
            gene_values = pd.read_csv(
                os.path.join(OUT_R_PATH, "Gene_Count", f"genes_cluster_{index}.csv"),
                sep="\t",
                engine="python",
            )["GENE"].values
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
        case 6:
            fig, ax = venn.venn6(labels, names=list_clusters, figsize=(10, 10))
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
    if len(list_clusters) < 2:
        return no_update
    CLUSTER_SELECTED_MULTI = list_clusters

    cluster_gene_list = {}
    all_gene_set = set()
    for index in list_clusters:
        gene_values = []
        if index == "ALL":
            gene_values = list(
                pd.read_csv(
                    os.path.join(OUT_R_PATH, "distribution_gene_cluster.csv"),
                    sep="\t",
                    engine="python",
                )["Gene"].unique()
            )
            cluster_gene_list["ALL"] = gene_values
        else:
            gene_values = list(
                pd.read_csv(
                    os.path.join(
                        OUT_R_PATH, "Gene_Count", f"genes_cluster_{index}.csv"
                    ),
                    sep="\t",
                    engine="python",
                )["GENE"].unique()
            )
            cluster_gene_list[index] = gene_values
        all_gene_set.update(gene_values)

    df = pd.DataFrame({"gene": list(all_gene_set)})
    for cluster, genes in cluster_gene_list.items():
        df[cluster] = df["gene"].apply(lambda x: "🟢" if x in genes else "🔴")

    df = df.sort_values(by=list_clusters, ascending=False)
    values = [df["gene"].values]
    c_name = ["Gene"]
    for i in list_clusters:
        c_name.append(f"Cluster {i}")
        values.append(df[i])

    fig = go.Figure(
        data=[
            go.Table(
                header={
                    "values": c_name,
                    "fill_color": "paleturquoise",
                    "align": "left",
                },
                cells={"values": values, "fill_color": "lavender", "align": "left"},
            )
        ]
    )
    return fig


# MULTI IMAGE BOX/PIE
def func_multi_plot(list_clusters, column_name):
    fig = None
    # Controllo dati
    df = DF_CLINICAL_DATA[DF_CLINICAL_DATA["cluster"].isin(list_clusters)]
    df.dropna(subset=[column_name], inplace=True)
    df["cluster_plot"] = df["cluster"].apply(lambda x: f"cluster_{x}")

    if column_name in NUMERIC_COLUMNS_CLINICAL:
        fig = tap.plot_stats(df, "cluster_plot", column_name)
    else:
        fig = make_subplots(
            1,
            len(list_clusters),
            subplot_titles=[f"cluster {e}" for e in list_clusters],
            specs=[[{"type": "domain"} for e in list_clusters]],
        )
        for i, index in enumerate(list_clusters):
            _cluster_values = df[df["cluster"] == index][column_name]
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
    Output("fig_multi_fig1", "figure"),
    [Input("dd-cluster-multi", "value"), Input("dd-box-1", "value")],
)
def update_multi_fig1(list_clusters, column_name):
    global CLUSTER_SELECTED_MULTI
    global BOX_FIG_SELECTED_1
    if len(list_clusters) < 2 or column_name is None:
        return no_update
    CLUSTER_SELECTED_MULTI = list_clusters
    BOX_FIG_SELECTED_1 = column_name
    return func_multi_plot(list_clusters, column_name)


# START
if __name__ == "__main__":
    print(f"{APP_NAME} ready on: http://127.0.0.1:8593")
    APP.run(debug=False, host="0.0.0.0", port=8593)
