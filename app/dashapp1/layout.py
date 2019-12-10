import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import pandas as pd
import metabolyze as met
from dash.dependencies import Input, Output
from functools import reduce
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import metabolyze as met
from functools import reduce
import warnings
import pandas as pd
import itertools
import scipy
import scipy.stats
import numpy as np
from functools import reduce
import re
import numpy 
import subprocess as sp
import os
import sys
import time
import dash
import dash_table
import pandas as pd
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input, State
import plotly.graph_objs as go
from sklearn.preprocessing import scale
import scipy.cluster.hierarchy as shc
import dash_table

def get_layout():

    result = met.Analysis('skeleton_output.tsv','Groups.csv')
    dropdown_comparison = [' vs. '.join(list(x)) for x in result.dme_comparisons()]
    #print("met_comparison",result.dme_comparisons())
    dropdown_list = []
    for dropdown,met_comparison in zip(dropdown_comparison,result.dme_comparisons()):
        drop_dict = {}
        drop_dict['label'] = dropdown
        met_list = list(met_comparison)
        met_list = ','.join("'{0}'".format(comparison) for comparison in met_list)
        drop_dict['value'] = met_list
        dropdown_list.append(drop_dict)


    logo = 'data:image/png;base64,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'



    search_bar = dbc.Row(
        [
            dbc.Col(dbc.Input(type="search", placeholder="Search")),
            dbc.Col(
                dbc.Button("Search", color="primary", className="ml-2"),
                width="auto",
            ),
        ],
        no_gutters=True,
        className="ml-auto flex-nowrap mt-3 mt-md-0",
        align="center",
    )

    navbar = dbc.Navbar(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src=logo, height="80")),
                        dbc.Col(dbc.NavbarBrand("Metabolive", className="ml-2")),
                    ],
                    align="center",
                    no_gutters=True,
                ),
                href="/login",

            ),
            dbc.NavItem(dbc.NavLink("Home", active=True, href="/",external_link=True)),
            dbc.NavItem(dbc.NavLink("Logout", active=True, href="/logout",external_link=True)),
            dbc.NavItem(dbc.NavLink("Contact Us", active=True, href="",external_link=True)),
            dbc.NavbarToggler(id="navbar-toggler"),

            dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
        ],

        color="dark",
        dark=True,
    )



    jumbotron = dbc.Jumbotron(
        [
            html.H1("Metabolive Beta", className="display-3"),
            html.P(
                "Beta testing for the interactive analysis of Metabolyze results."
            ),
        ]
    )

    tab1_content = dbc.Card(

    dbc.CardBody(
        [
            dbc.Button("Hybrid", size="lg",outline =True,color="primary"),
            dbc.Button("Untargeted", size="lg",outline =True,color="primary"),
            dbc.Button("Global", size="lg",outline =True,color="primary")
        ]
    ),
    className="mt-3",
)

    tab2_content = dbc.Card(
        dbc.CardBody(
            [
                dbc.Button("Run", outline=True,color="primary"),
            ]
        ),
        className="mt-3",
    )


    tabs = dbc.Tabs(
        [
            dbc.Tab(tab1_content, label="Metabolyze"),
            dbc.Tab(tab2_content, label="XIC Analysis")
        ]
    )
    #"'1e6 per mL','RPMI_PBS'"
    #[{'label': 'OLFR2 vs. GPI1', 'value': ['OLFR2', 'GPI1']}]

    button = html.Div(
        [
            dbc.Button("Submit", id="button"),
            html.Span(id="example-output"),
        ]
    )


    dropdown_processing = dcc.Dropdown(
            id='pipeline-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=[{'label':'Ungrid','value': 'Ungrid'}],
            value="Ungrid"
            #value=["test"]
        )



    dropdown_analysis = dcc.Dropdown(
            id='my-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=dropdown_list,
            #value=["test"]
        )
    table = html.Div(id='output-data-upload')

    hidden_div = html.Div(id='intermediate-value', style={'display': 'none'})
    volcano_graph = html.Div(id='volcano_graph')
    heatmap_graph = html.Div(id='heatmap_graph')
    pca_graph = html.Div(id='pca_graph')
    #graph = html.Div(id='graph')
    # def layout():
    layout = html.Div([navbar,jumbotron,tabs,dropdown_processing,dropdown_analysis,button,hidden_div,volcano_graph,heatmap_graph,pca_graph,table])
    return(layout)
    	# return(layout)

layout = (get_layout())
