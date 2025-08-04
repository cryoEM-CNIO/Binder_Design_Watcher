#!/usr/bin/env python3

'''
POSSIBLE PROBLEMS

- If you don't use the FastRelax cycle to improve the sequences, the names of the files change. In this case, you should remove the cycle_1 from the csv of Rosetta and CUTRE
This could have an easy fix

'''

'''
things to do: 

'''

import os
import pandas as pd

import dash
from dash import dash_table
import dash_bio as dashbio
from dash import html
import dash_bootstrap_components as dbc
from dash import Dash, dcc, Input, Output, State, callback
from dash.exceptions import PreventUpdate
import dash_bio.utils.ngl_parser as ngl_parser
import plotly.express as px
import subprocess
import re
from Bio import PDB
import argparse
import socket
import numpy as np
import time
import os.path
import math
from utils.hits_utils import *
from utils.generic_utils import *
from utils.plotting_utils import *




# Parsing port number and host
parser = argparse.ArgumentParser()
parser.add_argument("--port", "-p", help = "choose port to run the webapp")
parser.add_argument("--host", "-host", help = "choose host to run the webapp")
parser.add_argument("--debug", "-d", help = "launch app in debug mode")

args, unknown = parser.parse_known_args()

# Set localhost and 8051 as host and port by default
if not args.port: port_number = 8050
else: port_number = args.port
if not args.host: hostname = socket.gethostname()
else: hostname = args.host 
if not args.debug: debug_mode = False
else: debug_mode = True

# Messsssssinesssss
working_dir = os.getcwd()
output_dir = os.path.join(working_dir, 'output')
directories_list=get_working_directories(working_dir)
if not directories_list:
    directories_list=[working_dir]
designs_list=[]
print(directories_list)
#Initial lists for extraction dropdowns

initial_organisms=[
    "Arabidopsis thaliana",
    "Bacillus subtilis",
    "Caenorhabditis elegans",
    "Chlamydomonas reinhardtii",
    "Danio rerio",
    "Drosophila melanogaster",
    "Homo sapiens",
    "Mus musculus",
    "Nicotiana tabacum",
    "Pseudomonas putida",
    "Saccharomyces cerevisiae",
    "Escherichia coli general",
]

# Styles
bt_style = {"align-items": "center", "background-color": "#F2F3F4", "border": "2px solid #000",
            "box-sizing": "border-box", "color": "#000", "cursor": "pointer", "display": "inline-flex",
             'padding':'0.3em 1.2em', 'margin':'0.5em 0.3em 0.3em 0',
            "font-size": "0.9em", 'font-weight':'500', 'border-radius':'2em', 'text-align':'center',
            'transition':'all 0.2s'}
dropdown_style={'background-color':'#F2F3F4',
                "color": "#000", "cursor": "pointer",'width':'260px',
                 'margin':'0.5em 0.3em 0.3em 0',
                "font-size": "1.2em", 'font-weight':'500', 'text-align':'center',
                'transition':'all 0.2s'}

table_viewer_div_style={'display': 'flex', 'justifyContent': 'flex-start'}
table_div_style={'marginRight': '30px'}
table_style={'overflowX': 'auto','width': '100%', 'margin': '0 auto'}
table_cell_style={'minWidth': '150px', 'width': '150px', 'maxWidth': '150px', 'overflow': 'hidden', 'textOverflow': 'ellipsis', 'padding': '5px',}
title_style = {"margin-left": "15px", "margin-top": "15px", "margin-bottom": "0em", "color": "Black", "font-family" : "Helvetica", "font-size":"2.5em"}
box_style3 = {"font-size":"0.9em",'padding':'0.3em 1.2em','width':"18%", "margin-left": "0%","margin-right": "1%", "color": "black","font-family" : "Helvetica", 'vertical-align': 'center', "margin-bottom":"2px"}
extraction_box_style={'padding': '20px','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'height':'100px'}
extraction_box_cl_style={'padding': '20px','max-height':'240px','overflow-y':'auto','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'height':'290px'}
cool_button_style={'font-size': '24px','padding': '30px 60px','background': 'linear-gradient(135deg, #4CAF50, #81C784)','color': 'white','border': 'none','border-radius': '10px','box-shadow': '0px 4px 6px rgba(0, 0, 0, 0.1)','cursor': 'pointer','transition': 'transform 0.2s, box-shadow 0.2s','flex':'1','min-width': '300px', 'flex-basis': '350px','height': '290px'}
stop_button_style = {'font-size': '18px','padding': '30px 60px','background': 'linear-gradient(135deg, #FF0000, #FF6347)','color': 'white','border': 'none','border-radius': '10px','box-shadow': '0px 4px 6px rgba(0, 0, 0, 0.1)','cursor': 'pointer','transition': 'transform 0.2s, box-shadow 0.2s','flex': '1','min-width': '300px','flex-basis': '350px','height': '200px'}
extraction_cl_style={'margin-top':'50px','max-height':'750px','overflow-y':'auto','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'display':'flex', 'align-items':'flex-start'}
# Color definitions
color_points = 'cornflowerblue'

# Dash app initialization
# app = dash.Dash(__name__)

import dash_bootstrap_components as dbc

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)
app.title = "BinderFlow Monitor"
server = app.server

# Define layout
def serve_layout():
    return html.Div(
        className="app-wrapper",
        children=[
            # Sidebar
            html.Div(
                className="sidebar",
                children=[
                    html.H1("BinderFlow Monitor", className="app-title"),
                    # Removed: html.Div(id='row-count'),
                    dcc.Dropdown(
                        options=[
                            {
                                'label': (
                                    f"{os.path.basename(os.path.dirname(d))}/"
                                    f"{os.path.basename(d)}"
                                    if os.path.basename(os.path.dirname(d))
                                    else os.path.basename(d)
                                ),
                                'value': d
                            }
                            for d in directories_list
                        ],
                        value=directories_list[0],
                        placeholder='Select a folder',
                        id='directory-dropdown',
                        className='dropdown'
                    ),
                    html.Br(),
                    html.Button('STOP CAMPAIGN', id='stop-campaign', n_clicks=0, className="button button-danger"),
                    html.Br(),
                    html.Button('Make CSVs', id='execute-csv', n_clicks=0, className="button"),
                    html.Button('Filters & Axes', id='open-filters', n_clicks=0, className='button'),
                    dbc.Collapse(
                        html.Div([
                            html.H5('PAE Interaction'),
                            dcc.RangeSlider(id='pae_interaction_thres', min=0, max=30, value=[0,10], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('CUTRE'),
                            dcc.RangeSlider(id='CUTRE_thres', min=0, max=70, value=[0,10], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('pLDDT binder'),
                            dcc.RangeSlider(id='plddt_binder_thres', min=0, max=100, value=[80,100], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('dSASA'),
                            dcc.RangeSlider(id='dsasa_thres', min=0, max=10000, value=[1000,10000], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Shape Complementarity'),
                            dcc.RangeSlider(id='shape_complementarity_thres', min=0, max=1, value=[0.5,1], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Interface HBond'),
                            dcc.RangeSlider(id='interface_hbond_thres', min=0, max=15, value=[3,15], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Interface Unsaturated HBond'),
                            dcc.RangeSlider(id='interface_unsat_hbond_thres', min=0, max=15, value=[0,4], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Binder Surface Hydrophobicity'),
                            dcc.RangeSlider(id='binder_surf_hyd_thres', min=0, max=1, value=[0,0.35], tooltip={"placement":"bottom","always_visible":True}),
                            html.Hr(),
                            html.Div([
                                html.Span('X Value', className='axis-label'),
                                dcc.Dropdown( ['plddt_binder','pae_interaction','CUTRE','dG','dSASA','Shape_complementarity','Packstat','dG_SASA_ratio','length','SAP','binder_int_hyd','binder_surf_hyd','interface_hbonds','interface_unsat_hbonds','ipSAE','RMSD'], 'pae_interaction', id='xaxis_value', className='dropdown')
                            ], className='axis-control'),
                            html.Div([
                                html.Span('Y Value', className='axis-label'),
                                dcc.Dropdown( ['plddt_binder','pae_interaction','CUTRE','dG','dSASA','Shape_complementarity','Packstat','dG_SASA_ratio','length','SAP','binder_int_hyd','binder_surf_hyd','interface_hbonds','interface_unsat_hbonds','ipSAE','RMSD'], 'plddt_binder', id='yaxis_value', className='dropdown')
                            ], className='axis-control'),
                        ], style={'padding':'10px'}),
                        id='filters-collapse',
                        is_open=False
                    ),
                    html.Br(),
                ]
            ),
            # Main content
            html.Div(
                className="main-content",
                children=[
                    html.Div([
                        dcc.Tabs([
                            dcc.Tab(label='Live Watcher', children=[
                                # Metrics cards row
                                html.Div(
                                    children=[
                                        dbc.Card(
                                            [
                                                html.Div(id='row-count', className='metric-value'),
                                                html.Div("Finished Designs", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                        dbc.Card(
                                            [
                                                html.Div(id='hit-count', className='metric-value'),
                                                html.Div("Hits", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                        dbc.Card(
                                            [
                                                html.Div(id='hit-efficiency', className='metric-value'),
                                                html.Div("Hit Efficiency", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                    ],
                                    className="metrics-row"
                                ),
                                html.Div([
                                    html.Div([
                                        html.Div([
                                            # Removed filter/axis toggles and panels; controls now in Offcanvas
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader("Scatter Plot"),
                                                    dbc.CardBody(
                                                        dcc.Graph(
                                                            id='scatter-plot',
                                                            className='graph-container',
                                                            config={'toImageButtonOptions': {'format': 'svg', 'filename': 'scatter_plot', 'height': 600, 'width': 800, 'scale': 1}, 'displaylogo': False}
                                                        )
                                                    ),
                                                ],
                                                className="graph-card"
                                            ),
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader("Radar Plot"),
                                                    dbc.CardBody(
                                                        dcc.Graph(
                                                            id='radar-plot',
                                                            className='graph-container',
                                                            config={'toImageButtonOptions': {'format': 'svg', 'filename': 'radar_plot', 'height': 600, 'width': 800, 'scale': 1}, 'displaylogo': False}
                                                        )
                                                    ),
                                                ],
                                                className="graph-card"
                                            ),
                                        ], className='plot-and-controls'),
                                    ], className='plot-and-controls'),
                                ]),
                                html.Div(
                                    [
                                        html.Div(
                                            [
                                                html.H4('Original Input Path (For Partial Diffusion comparison)'),
                                                dcc.Input(
                                                    id='input_pdb_path',
                                                    type='text',
                                                    placeholder='Input PDB file',
                                                    className="input-box"
                                                )
                                            ],
                                            style={'width':'100%', 'padding':'10px'}
                                        )
                                    ],
                                    style={'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'}
                                ),
                            ]),
                            dcc.Tab(label='Pipeline Tracking', children=[
                                html.Div(id='job-status-counts', style={'margin': '20px 0'}),
                                dbc.ListGroup(id='job-status-list', flush=True, className='mb-4'),
                                dcc.Interval(id='interval-component', interval=60000, n_intervals=0),
                                dcc.Store('filtered_df')
                            ]),
                        dcc.Tab(label='Extraction', children=[
                            # Hit Viewer & Selection
                            dbc.Row([
                                dbc.Col(
                                    dbc.Card([
                                        dbc.CardHeader("Hit PDB preview", className='card-header-primary'),
                                        dbc.CardBody([
                                            dcc.Dropdown(
                                                options=[],
                                                value=None,
                                                placeholder='Select a hit',
                                                id='extractions_molecule_dropdown',
                                                className='dropdown'
                                            ),
                                            dashbio.NglMoleculeViewer(id="extractions_molecule")
                                        ])
                                    ], className='mb-4'),
                                    width=4
                                ),
                                dbc.Col(
                                    dbc.Card([
                                        dbc.CardHeader("Selection of hits for extraction", className='card-header-primary'),
                                        dbc.CardBody([
                                            dcc.Checklist(
                                                options=[],
                                                value=[],
                                                id='extraction-selection',
                                                style={'display':'flex','flexDirection':'column','gap':'5px'}
                                            )
                                        ])
                                    ], className='mb-4'),
                                    width=8
                                )
                            ]),
                            # Extraction Options Card
                            dbc.Card([
                                dbc.CardHeader("Extraction Options", className='card-header-primary'),
                                dbc.CardBody([
                                    dbc.Row([
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("Source File"),
                                                dbc.CardBody(
                                                    dcc.RadioItems(
                                                        ['PMPNN', 'AF2'],
                                                        'AF2',
                                                        id='extract-type',
                                                        inputStyle={'margin-right':'5px'}
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        ),
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("Generate DNA Sequence?"),
                                                dbc.CardBody(
                                                    dcc.RadioItems(
                                                        ['Yes', 'No'],
                                                        'Yes',
                                                        id='DNA-seq',
                                                        inputStyle={'margin-right':'5px'}
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        )
                                    ]),
                                    dbc.Row([
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("Add Initial Methionine"),
                                                dbc.CardBody(
                                                    dcc.RadioItems(
                                                        ['True', 'False'],
                                                        'False',
                                                        id='add_met',
                                                        inputStyle={'margin-right':'5px'}
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        ),
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("Organism"),
                                                dbc.CardBody(
                                                    dcc.Dropdown(
                                                        initial_organisms,
                                                        'Escherichia coli general',
                                                        id='organism',
                                                        className='dropdown'
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        )
                                    ]),
                                    dbc.Row([
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("3' Overhang"),
                                                dbc.CardBody(
                                                    dcc.Input(
                                                        placeholder='Sequence',
                                                        value='',
                                                        id='three_prime_overhang',
                                                        className='input-box'
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        ),
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("5' Overhang"),
                                                dbc.CardBody(
                                                    dcc.Input(
                                                        placeholder='Sequence',
                                                        value='',
                                                        id='five_prime_overhang',
                                                        className='input-box'
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        )
                                    ]),
                                    dbc.Row([
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("Random Sequence Length"),
                                                dbc.CardBody(
                                                    dcc.Input(
                                                        type='number',
                                                        value=300,
                                                        id='random_sequence',
                                                        className='input-box'
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        ),
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("GC % of Random Sequence"),
                                                dbc.CardBody(
                                                    dcc.Input(
                                                        type='number',
                                                        value=50,
                                                        id='GC_content',
                                                        className='input-box'
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=6
                                        )
                                    ]),
                                    dbc.Row([
                                        dbc.Col(
                                            dbc.Card([
                                                dbc.CardHeader("Check Restriction Sites"),
                                                dbc.CardBody(
                                                    dcc.Checklist(
                                                        options=[{'label': enz, 'value': enz} for enz in [
                                                            'EcoRI','BamHI','HindIII','NotI','XhoI','PstI','SacI','KpnI',
                                                            'SmaI','XbaI','SpeI','NcoI','SalI','ApaI','HaeIII','AluI',
                                                            'TaqI','BglII','ClaI','MluI','BsaI'
                                                        ]],
                                                        value=[],
                                                        id='enzyme',
                                                        inputStyle={'margin-right':'5px'},
                                                        className='enzyme-grid'
                                                    )
                                                )
                                            ], className='mb-3'),
                                            width=12
                                        )
                                    ]),
                                    html.Div(
                                        dbc.Button('Extract Hits', id='execute-hits', n_clicks=0, className='button'),
                                        style={'textAlign':'center', 'margin-top':'10px'}
                                    )
                                ])
                            ], className='mb-4')
                        ]),
                        ])
                    ])
                ]
            ),
            # Offcanvas removed
        ]
    )

app.layout = serve_layout


## PDB viewer
@callback(
    Output("extractions_molecule", 'data'),
    Output("extractions_molecule", "molStyles"),
    Input("extractions_molecule_dropdown", "value"),
    Input('directory-dropdown', 'value')
)

def return_molecule(value,directory):
    working_dir=f'{directory}/output/'
    if (value is None) or (value == ''):
        raise PreventUpdate
    else:
        data_path, filename = get_design_file_path_and_name(value, working_dir)
        try: data_path = data_path+'/'
        except: data_path = ''
        molstyles_dict = {
            "representations": ["ribbon"],
            "chosenAtomsColor": "white",
            "chosenAtomsRadius": 1,
            "molSpacingXaxis": 100,
        }
    data_list = [
        # chain A in blue, reset_view=True to recenter on the first model
        ngl_parser.get_data(
            data_path=data_path,
            pdb_id=f"{filename}.A",    # tells NGL to load only chain A
            color='#edb081',
            reset_view=True,
            local=True
        ),
        # chain B in green, reset_view=False to keep the same camera
        ngl_parser.get_data(
            data_path=data_path,
            pdb_id=f"{filename}.B",    # load only chain B
            color='#4b2362',
            reset_view=False,
            local=True
        )
    ]

    return data_list, molstyles_dict

# Callback to update graphs and table
@callback(
    Output('scatter-plot', 'figure'),
    Output('row-count', 'children'),
    Output('hit-count', 'children'),
    Output('hit-efficiency', 'children'),
    Output('job-status-counts', 'children'),
    Output('extractions_molecule_dropdown', 'options'),
    Output('filtered_df', 'data'),
    Output('extraction-selection', 'options'),
    Output('extraction-selection', 'value'),
    [
        Input('directory-dropdown', 'value'),
        Input('interval-component', 'n_intervals'),
        Input('stop-campaign', 'n_clicks'),
        Input('xaxis_value', 'value'),
        Input('yaxis_value','value'),
        Input('input_pdb_path', 'value'),
        Input('pae_interaction_thres', 'value'),
        Input('CUTRE_thres', 'value'),
        Input('plddt_binder_thres', 'value'),
        Input('dsasa_thres', 'value'),
        Input('shape_complementarity_thres', 'value'),
        Input('interface_hbond_thres', 'value'),
        Input('interface_unsat_hbond_thres', 'value'),
        Input('binder_surf_hyd_thres', 'value'),
    ]
)

def update_graph(working_dir, n, n_clicks_stop, xaxis_value, yaxis_value, input_pdb_path, pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    directory = f'{working_dir}/output/'
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(working_dir, input_pdb_path)
    filtered_df = filtering_df(merged_df, pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres)
    status_df = pd.DataFrame()
    status_df = track_job_status(directory)
    status_df = status_df
    if not status_df.empty:
        status_df_records = status_df.to_dict('records')
    else:
        status_df_records = []

    # Make scatter plot
    scatter_plot, row_count_text = update_scatter_plot(working_dir, merged_df, filtered_df, xaxis_value, yaxis_value, input_pdb_path)

    ctx = dash.callback_context

    # STOP CAMPAIGN
    if ctx.triggered and ctx.triggered[0]['prop_id'] == 'stop-campaign.n_clicks':
        command = f'touch {working_dir}/campaign_done'
        subprocess.run(command, shell=True)

    # Job status counts
    status_counts = status_df['status'].value_counts().to_dict()
    status_counts_formatted = ", ".join([f"{status}: {count}" for status, count in status_counts.items()])
    job_status_counts_text = f"Job Status Counts: {status_counts_formatted}"
    status_counts = 0
    job_status_counts_text = f"Job Status Counts: 0"

    # Dropdown HITS
    dropdown_options = get_hit_names(filtered_df, xaxis_value)

    # jasonified df for communication between callbacks
    jasonified_df = filtered_df.to_json(date_format='iso', orient='split')

    # Metrics for cards
    finished_models = len(merged_df)
    hit_count = len(dropdown_options)
    if finished_models > 0:
        hit_efficiency = f"{(hit_count/finished_models*100):.1f}%"
    else:
        hit_efficiency = "N/A"

    return (
        scatter_plot,
        row_count_text,
        hit_count,
        hit_efficiency,
        job_status_counts_text,
        dropdown_options,
        jasonified_df,
        dropdown_options,
        dropdown_options
    )
# New callback for job-status-list
@callback(
    Output('job-status-list', 'children'),
    Input('interval-component', 'n_intervals'),
    State('directory-dropdown', 'value')
)
def update_job_list(n, working_dir):
    df = track_job_status(f'{working_dir}/output/')
    # Extract run number and sort descending
    df_sorted = df.copy()
    df_sorted['run_number'] = df_sorted['job'].str.extract(r'(\d+)', expand=False).astype(int)
    df_sorted = df_sorted.sort_values('run_number', ascending=False)
    items = []
    for _, row in df_sorted.iterrows():
        status = row['status'].upper()
        color = {
            'RUNNING': 'warning',
            'SUCCESS': 'success',
            'FAILED': 'danger',
            'WAITING': 'info'
        }.get(status, 'secondary')
        items.append(
            dbc.ListGroupItem(
                [
                    html.Span(row['job'], className='me-2', style={'fontWeight': '600'}),
                    dbc.Badge(status, color=color, pill=True, className='ms-auto')
                ],
                className='d-flex align-items-center'
            )
        )
    return items

@callback(
    Output('radar-plot', 'figure'),
    [
    Input('scatter-plot', 'clickData'),
    Input('interval-component', 'n_intervals'),
    Input('directory-dropdown', 'value'),
    Input('input_pdb_path', 'value'),
    Input('pae_interaction_thres', 'value'),
    Input('CUTRE_thres', 'value'),
    Input('plddt_binder_thres', 'value'),
    Input('dsasa_thres', 'value'),
    Input('shape_complementarity_thres', 'value'),
    Input('interface_hbond_thres', 'value'),
    Input('interface_unsat_hbond_thres', 'value'),
    Input('binder_surf_hyd_thres', 'value'),
    ])

#This prints the radar plot
def update_radar_plot(design_to_plot, n, directory, input_pdb_path, pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    #This is meant to store the original designs for the following plotting; is a little bit cutre 
    update_designs_list(designs_list, design_to_plot)
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(directory, input_pdb_path)
    #Most of the heavy work is carry in the utils file, go there for further info
    radar_figure=radar_plot(designs_list, merged_df,pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres )
    # Override radar area colors with transparency
    fill_colors = ['rgba(134,48,113,0.5)', 'rgba(237,176,129,0.5)']
    line_colors = ['#863071', '#edb081']
    for i, trace in enumerate(radar_figure.data):
        if hasattr(trace, 'fillcolor'):
            trace.fillcolor = fill_colors[i % len(fill_colors)]
            trace.line.color = line_colors[i % len(line_colors)]
        elif hasattr(trace, 'marker'):
            trace.marker.color = line_colors[i % len(line_colors)]
    return radar_figure

@callback(
    [
        Input('directory-dropdown', 'value'),
        Input('interval-component', 'n_intervals'),
        Input('execute-hits', 'n_clicks'),
        Input('extract-type','value'),
        Input('DNA-seq', 'value'),
        Input('add_met', 'value'),
        Input('organism', 'value'),
        Input('three_prime_overhang', 'value'),
        Input('five_prime_overhang', 'value'),
        Input('extraction-selection', 'value'),
        Input('random_sequence', 'value'),
        Input('GC_content', 'value'),
        Input('enzyme', 'value')
    ]
)

def extract_hits(working_dir, n, clicks, extract_type,DNA_seq, met, organism, three_prime, five_prime, extraction_list, length, GC,enzyme):
    
    #Create the hits folder
    hits_folder = os.path.join(working_dir, 'hits')
    fastas_folder = os.path.join(hits_folder, 'fastas')
    dnaseq_folder=os.path.join(hits_folder, 'dna_seqs')
    if not os.path.exists(hits_folder):
        os.makedirs(hits_folder)
        os.makedirs(fastas_folder)
        os.makedirs(dnaseq_folder)
    
    ctx = dash.callback_context

    if ctx.triggered and ctx.triggered[0]['prop_id'] == 'execute-hits.n_clicks':
        if extraction_list:
            for description in extraction_list:
                print('###############################\nEXTRACTING HIT\n###############################\n' + description)
                run_number = re.search(r'run_\d+',description).group()
                design_number = math.floor(int(re.search(r'run_\d+_design_(\d+).*',description).group(1))/10)
                ########## For the v2 tests microrun, different structure ######################################################
                ########## If we go forward with the V3, marked lines must be remove before deployment ######################### 
                if design_number != 0:
                    pmpnn_file = f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_out.silent'
                    af2_file_sol=f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_sol_out_af2.silent'
                    af2_file=f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_out_af2.silent'
                else:
                    if os.path.isfile(f'{working_dir}/output/{run_number}/{run_number}_input_out.silent'):
                        pmpnn_file = f'{working_dir}/output/{run_number}/{run_number}_input_out.silent'
                        af2_file_sol=f'{working_dir}/output/{run_number}/{run_number}_input_sol_out_af2.silent'
                        af2_file=f'{working_dir}/output/{run_number}/{run_number}_input_out_af2.silent'
                    else:
                        design_number=int(re.search(r'run_\d+_design_(\d+).*',description).group(1))
                        pmpnn_file = f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_out.silent'
                        af2_file_sol=f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_sol_out_af2.silent'
                        af2_file=f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_out_af2.silent'
                if extract_type == 'PMPNN':
                        command = f'silentextractspecific {pmpnn_file}' + description[:-8] + ' > extraction.log'

                else:

                    if os.path.isfile(af2_file_sol):
                        command = f'silentextractspecific {af2_file_sol} '+ description + ' > extraction.log'  
                        subprocess.run(command, cwd=hits_folder, shell=True)
                    elif os.path.isfile(af2_file):
                        command = f'silentextractspecific {af2_file} ' + description + ' > extraction.log'  
                        subprocess.run(command, cwd=hits_folder, shell=True)
                    else:
                        print('The AF2 file required is not being found...')
                        continue
                #Get the fastas of the hits
                fasta_file=extract_fasta_seq(os.path.join(hits_folder, description))
                if DNA_seq == 'Yes':
                    extract_dna_seq(input=fasta_file,output=dnaseq_folder, organism=organism, met=met,overhang_3=three_prime, overhang_5=five_prime, length=length, GC=GC,enzyme=enzyme )
                #Record scoring data in the pdb
                add_stats_to_pdb(description, working_dir)
                print('File extracted')
            if DNA_seq == 'Yes':
                generate_order_csv(extraction_list,hits_folder)
                create_log_extraction(hits_folder, extraction_list, met, organism, three_prime, five_prime, length, GC)




# Callback to toggle the Collapse for filters and axes in the sidebar
@callback(
    Output('filters-collapse', 'is_open'),
    Input('open-filters', 'n_clicks'),
    State('filters-collapse', 'is_open')
)
def toggle_filters_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
# Run the app
if __name__ == '__main__':
    app.run(debug=False, dev_tools_hot_reload = False, use_reloader=True,
                   host=hostname, port=port_number)
