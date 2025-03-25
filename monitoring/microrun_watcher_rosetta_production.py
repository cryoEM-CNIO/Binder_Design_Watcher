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
from dash import Dash, dcc, html, Input, Output, callback
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

# three_prime_overhangs=[
#     'EGFP-10xHis', 
#     '6xHis-tag', 
#     'GST-tag', 
#     'FLAG-tag', 
#     'MBP-tag', 
#     'T7-tag', 
#     'Twin-Strep-tag',
#     'C-terminal_linker'
# ]

# five_prime_overhangs=[
#     'EGFP-10xHis', 
#     '6xHis-tag', 
#     'GST-tag', 
#     'FLAG-tag', 
#     'MBP-tag', 
#     'T7-tag', 
#     'Twin-Strep-tag', 
#     'N-terminal_linker', 
#     'Gateway'
# ]

# Styles
bt_style = {"align-items": "center", "background-color": "#F2F3F4", "border": "2px solid #000",
            "box-sizing": "border-box", "color": "#000", "cursor": "pointer", "display": "inline-flex",
            "font-family": "Helvetica", 'padding':'0.3em 1.2em', 'margin':'0.5em 0.3em 0.3em 0',
            "font-size": "0.9em", 'font-weight':'500', 'border-radius':'2em', 'text-align':'center',
            'transition':'all 0.2s'}
dropdown_style={'background-color':'#F2F3F4',
                "color": "#000", "cursor": "pointer",'width':'260px',
                "font-family": "Helvetica", 'margin':'0.5em 0.3em 0.3em 0',
                "font-size": "1.2em", 'font-weight':'500', 'text-align':'center',
                'transition':'all 0.2s'}

table_viewer_div_style={'display': 'flex', 'justifyContent': 'flex-start'}
table_div_style={'marginRight': '30px'}
table_style={'overflowX': 'auto','width': '100%', 'margin': '0 auto'}
table_cell_style={'minWidth': '150px', 'width': '150px', 'maxWidth': '150px', 'overflow': 'hidden', 'textOverflow': 'ellipsis', 'padding': '5px',}
title_style = {"margin-left": "15px", "margin-top": "15px", "margin-bottom": "0em", "color": "Black", "font-family" : "Helvetica", "font-size":"2.5em"}
box_style3 = {"font-size":"0.9em",'padding':'0.3em 1.2em','width':"18%", "margin-left": "0%","margin-right": "1%", "color": "black","font-family" : "Helvetica", 'vertical-align': 'center', "margin-bottom":"2px"}
extraction_box_style={'padding': '20px','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px'}
cool_button_style={'font-size': '18px','padding': '30px 60px','background': 'linear-gradient(135deg, #4CAF50, #81C784)','color': 'white','border': 'none','border-radius': '10px','box-shadow': '0px 4px 6px rgba(0, 0, 0, 0.1)','cursor': 'pointer','transition': 'transform 0.2s, box-shadow 0.2s','flex':'1','min-width': '300px', 'flex-basis': '350px','height': '200px'}
# Color definitions
color_points = 'cornflowerblue'

# Dash app initialization
# app = dash.Dash(__name__)

app = dash.Dash(
    __name__,
    assets_folder=output_dir,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)
app.title = "RFD Microruns watcher"
server = app.server

# Define layout
def serve_layout(): 
    return html.Div([

    #Title
    html.H1("RFD Microruns watcher", style=title_style),
    html.Div(id='row-count'),  
    dcc.Dropdown(options=directories_list, value=directories_list[0], placeholder='Select a folder', id='directory-dropdown'),

    #Tabs
    html.Div([
        dcc.Tabs([
            dcc.Tab(label='Live Watcher', children=[
                html.Button('Make CSVs', id='execute-csv', n_clicks=0, style=bt_style),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H5('X Value'),
                                dcc.Dropdown(
                                    ['plddt_binder','pae_interaction','CUTRE','dG','dSASA','Shape_complementarity','Packstat','dG_SASA_ratio','length','SAP','binder_int_hyd','binder_surf_hyd','interface_hbonds','interface_unsat_hbonds'],
                                    'pae_interaction',
                                    id='xaxis_value',
                                    style=dropdown_style
                                ),
                            ],
                            style={'display': 'inline-block', 'verticalAlign': 'top'}
                        ),
                        html.Div(
                            [
                                html.H5('Y Value'),
                                dcc.Dropdown(
                                    ['plddt_binder','pae_interaction','CUTRE','dG','dSASA','Shape_complementarity','Packstat','dG_SASA_ratio','length','SAP','binder_int_hyd','binder_surf_hyd','interface_hbonds','interface_unsat_hbonds'],
                                    'plddt_binder',
                                    id='yaxis_value',
                                    style=dropdown_style
                                ),
                            ],
                            style={'display': 'inline-block', 'verticalAlign': 'top'}
                        ),
                    ],
                    style={'display': 'flex', 'alignItems': 'center'}
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                dcc.Graph(id='scatter-plot', style={'width': '80vh', 'height': '70vh'}),
                                html.Div([

                                ])
                            ]
                        ),
                        html.Div(
                            [
                                html.H5('PAE Interaction thresholds'),
                                dcc.RangeSlider(
                                    id='pae_interaction_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=30,     # Define the maximum value of the slider
                                    value=[0, 10], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.H5('CUTRE thresholds'),
                                dcc.RangeSlider(
                                    id='CUTRE_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=70,     # Define the maximum value of the slider
                                    value=[0, 10], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.H5('pLDDT binder thresholds'),
                                dcc.RangeSlider(
                                    id='plddt_binder_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=100,     # Define the maximum value of the slider
                                    value=[80, 100], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.H5('dSASA thresholds'),
                                dcc.RangeSlider(
                                    id='dsasa_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=10000,     # Define the maximum value of the slider
                                    value=[1000, 10000], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.H5('Shape Complementarity thresholds'),
                                dcc.RangeSlider(
                                    id='shape_complementarity_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=1,     # Define the maximum value of the slider
                                    value=[0.5, 1], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.H5('Interface HBond thresholds'),
                                dcc.RangeSlider(
                                    id='interface_hbond_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=15,     # Define the maximum value of the slider
                                    value=[3, 15], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.H5('Interface Unsaturated HBond thresholds'),
                                dcc.RangeSlider(
                                    id='interface_unsat_hbond_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=15,     # Define the maximum value of the slider
                                    value=[0, 4], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.H5('Binder Surface Hydrophobicity thresholds'),
                                dcc.RangeSlider(
                                    id='binder_surf_hyd_thres',
                                    min=0,       # Define the minimum value of the slider
                                    max=1,     # Define the maximum value of the slider
                                    value=[0, 0.35], # Initial range
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                            ],
                            style={'width': '30vh', 'padding': '0 10px'}  # Adjust width for alignment

                        ),
                    
                        html.Div(
                            [
                                dcc.Graph(id='radar-plot', style={'width': '80vh', 'height': '70vh'}),
            
                            ]
                        ),
                    ],
                    style={'display': 'flex','align-items': 'center' }
                ),
                html.Div(
                    [
                        html.H4('Original Input Path (For Partial Diffusion comparison)'),
                        dcc.Input(id='input_pdb_path', type='text', placeholder='Input PDB file', style=box_style3)
                    ],
                    style={'width':'100%', 'padding':'10px'}
                ),
                html.Div(id='job-status-counts', style={'margin': '20px 0'}),
                html.Div([
                    html.Div([
                        dash_table.DataTable(
                            id='job-status-table',
                            columns=[
                                {"name": "Sort Key", "id": "sort_key"},
                                {"name": "Job Name", "id": "job", "deletable": False, "selectable": False},
                                {"name": "Run Status", "id": "status", "deletable": False, "selectable": False},
                            ],
                            data=[],
                            sort_action="native",
                            filter_action="native",
                            page_size=20,  # Set number of rows per page
                            style_cell=table_cell_style,
                            style_table=table_style,
                            )],style=table_div_style),
                    ],style=table_viewer_div_style),
                dcc.Interval(
                id='interval-component',
                interval=60000,  # in milliseconds, update every 20 seconds
                n_intervals=0
                ),
                dcc.Store('filtered_df')
            ]),
            dcc.Tab(label='Extraction', children=[
                html.Div([
                    html.Div([
                        dcc.Dropdown(options='', value='', placeholder='Select a hit', id='extractions_molecule_dropdown'),
                        dashbio.NglMoleculeViewer(id="extractions_molecule")
                            ], style={'width':'50%', 'padding':'10px'}
                        ),
                    html.Div([
                        dcc.Checklist(options='', value='', id='extraction-selection' )
                        ], style={'width':'50%', 'padding':'10px'})
                    ], style={'display': 'flex','align-items': 'center' },
                ),
                html.Div([
                    html.H3("Extraction options"),
                    html.Div([
                            html.Div([
                                html.H4('File to extract from:', style={'margin-bottom': '5px'}),
                                dcc.RadioItems(
                                    ['PMPNN', 'AF2'], 
                                    'AF2', 
                                    id='extract-type', 
                                    style={'margin-bottom': '20px', 'display': 'block'}
                                ),
                            ], style=extraction_box_style), 

                            html.Div([
                                html.H4("Add initial Methionine:", style={'margin-bottom': '5px'}),
                                dcc.RadioItems(
                                    ['True', 'False'], 
                                    'False', 
                                    id='add_met', 
                                    style={'display': 'block'}
                                ),
                            ], style=extraction_box_style),

                            html.Div([
                                html.H4("Organism"),
                                dcc.Dropdown(
                                    initial_organisms,
                                    'Escherichia coli general',
                                    id='organism',
                                    style={'display':'block'})
                                ], style=extraction_box_style),
                            html.Div([       
                                html.H4("3\' overhang"),
                                dcc.Input(
                                    placeholder='Input your sequence',
                                    value='',
                                    id='three_prime_overhang',
                                    style={'display':'block'})
                                ], style=extraction_box_style),
                            html.Div([
                                html.H4("5\' overhang"),
                                dcc.Input(
                                    placeholder='Input your sequence',
                                    value='',
                                    id='five_prime_overhang',
                                    style={'display':'block'})
                                    ], style=extraction_box_style),
                            html.Div([
                                html.H4("Reach a length of ..."),
                                dcc.Input(
                                    value=300,
                                    type="number",
                                    id='random_sequence',
                                    style={'display':'block'})
                                    ], style=extraction_box_style),     
                            html.Div([
                                html.H4("GC % of the random sequence"),
                                dcc.Input(
                                    value=50,
                                    id='GC_content',
                                    type="number",
                                    style={'display':'block'})
                                    ]),                         
                            html.Div([
                                html.Button(
                                    'Extract hits',
                                    id='execute-hits',
                                    n_clicks=0,
                                    style=cool_button_style,
                                    )
                                ]),
                        html.Div(id='hidden-output', style={'display': 'none'})
                        ], style={
                            'display': 'flex', 
                            'flex-direction': 'row', 
                            'gap': '10px', 
                            'width': '60vh', 
                            'padding': '0 10px'
                            }),

                    ])
                ])                   
            ])
        ])
    ])
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
            "representations": ["cartoon"],
            "chosenAtomsColor": "white",
            "chosenAtomsRadius": 1,
            "molSpacingXaxis": 100,
        }
        data_list = [ngl_parser.get_data(data_path=data_path, pdb_id=filename, color='red',reset_view=True, local=True)]
    return data_list, molstyles_dict

# Callback to update graphs and table
@callback(
    [Output('scatter-plot', 'figure'),
     Output('row-count', 'children'),
     Output('job-status-table', 'data'),
     Output('job-status-counts', 'children'),
     Output('extractions_molecule_dropdown', 'options'),
     Output('filtered_df', 'data'),
     Output('extraction-selection', 'options'),
     Output('extraction-selection', 'value')],
    [Input('directory-dropdown', 'value'),
     Input('interval-component', 'n_intervals'),
     Input('execute-csv', 'n_clicks'),
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
     ])

def update_graph( working_dir,n, n_clicks_csv,xaxis_value,yaxis_value, input_pdb_path, pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
        
    directory=f'{working_dir}/output/'
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(working_dir, input_pdb_path)
    filtered_df=filtering_df(merged_df, pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres )
    status_df = pd.DataFrame()
    status_df = track_job_status(directory)
    status_df=status_df
    if not status_df.empty:
        status_df_records = status_df.to_dict('records')
    else:
        status_df_records = []

    ctx = dash.callback_context
    
    #Generate CSVs
    if ctx.triggered[0]['prop_id'] == 'execute-csv.n_clicks':
        if not merged_df.empty:
            merged_df.iloc[:, 1:].to_csv(f'{working_dir}/all_models.csv', index=False)
            if not filtered_df.iloc[:, 1:].empty:
                filtered_df.to_csv(f'{working_dir}/hits.csv', index=False)
        scatter_plot, row_count_text = update_scatter_plot(working_dir,merged_df,filtered_df, xaxis_value, yaxis_value,input_pdb_path )
    else:
        scatter_plot, row_count_text = update_scatter_plot(working_dir, merged_df,filtered_df,xaxis_value, yaxis_value, input_pdb_path )
    # Job status counts
    status_counts = status_df['status'].value_counts().to_dict()
    status_counts_formatted = ", ".join([f"{status}: {count}" for status, count in status_counts.items()])
    job_status_counts_text = f"Job Status Counts: {status_counts_formatted}"
    status_counts=0
    job_status_counts_text = f"Job Status Counts: 0"


    # Dropdown HITS
    
    dropdown_options = get_hit_names(filtered_df,xaxis_value)

    #jasonified df for communication between callbacks
    jasonified_df=filtered_df.to_json(date_format='iso', orient='split')


    return scatter_plot, row_count_text, status_df_records, job_status_counts_text,dropdown_options, jasonified_df, dropdown_options, dropdown_options

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
    ]
)

#This prints the radar plot
def update_radar_plot(design_to_plot, n, directory, input_pdb_path, pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    #This is meant to store the original designs for the following plotting; is a little bit cutre 
    update_designs_list(designs_list, design_to_plot)
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(directory, input_pdb_path)
    #Most of the heavy work is carry in the utils file, go there for further info
    radar_figure=radar_plot(designs_list, merged_df,pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres )
    return radar_figure

@callback(
    Output('hidden-output', 'children'),
    [
        Input('directory-dropdown', 'value'),
        Input('interval-component', 'n_intervals'),
        Input('execute-hits', 'n_clicks'),
        Input('extract-type','value'),
        Input('add_met', 'value'),
        Input('organism', 'value'),
        Input('three_prime_overhang', 'value'),
        Input('five_prime_overhang', 'value'),
        Input('extraction-selection', 'value'),
        Input('random_sequence', 'value'),
        Input('GC_content', 'value')
    ]
)

def extract_hits(working_dir, n, clicks, extract_type, met, organism, three_prime, five_prime, extraction_list, length, GC):
    
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
                if extract_type == 'PMPNN':
                    command = f'silentextractspecific {working_dir}/output/'+ run_number + '/' + run_number + '_input_out.silent ' + description[:-8] + ' > extraction.log'
                else:
                    command = f'silentextractspecific {working_dir}/output/'+ run_number + '/' + run_number + '_input_sol_out_af2.silent ' + description + ' > extraction.log'  
                subprocess.run(command, cwd=hits_folder, shell=True)
                #Record scoring data in the pdb
                add_stats_to_pdb(description, working_dir)
                
                #Get the fastas of the hits
                fasta_file=extract_fasta_seq(os.path.join(hits_folder, description))
                extract_dna_seq(input=fasta_file,output=dnaseq_folder, organism=organism, met=met,overhang_3=three_prime, overhang_5=five_prime, length=length, GC=GC )
            generate_order_csv(extraction_list)
            create_log_extraction(hits_folder, extraction_list, met, organism, three_prime, five_prime, length, GC)
    return None

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, dev_tools_hot_reload = False, use_reloader=True,
                   host=hostname, port=port_number)
