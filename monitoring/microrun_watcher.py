#!/usr/bin/env python3

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

# Parsing port number and host
parser = argparse.ArgumentParser()
parser.add_argument("--port", "-p", help = "choose port to run the webapp")
parser.add_argument("--host", "-host", help = "choose host to run the webapp")
parser.add_argument("--debug", "-d", help = "launch app in debug mode")
parser.add_argument("--pae_thres", "-pae", type=float, help = "threshold for pae_interaction for extracting hits")
parser.add_argument("--plddt_thres", "-plddt", type=float, help = "threshold for plddt_binder for extracting hits")

args, unknown = parser.parse_known_args()

# Set localhost and 8051 as host and port by default
if not args.port: port_number = 8050
else: port_number = args.port
if not args.host: hostname = socket.gethostname()
else: hostname = args.host 
if not args.debug: debug_mode = False
else: debug_mode = args.host
if not args.pae_thres: pae_thres = 10
else: pae_thres = args.pae_thres
if not args.plddt_thres: plddt_thres = 80
else: plddt_thres = args.plddt_thres

# Messsssssinesssss
working_dir = os.getcwd()
output_dir = os.path.join(working_dir, 'output')
hits_folder = os.path.join(os.getcwd(), 'hits')
if not os.path.exists(hits_folder):
    os.makedirs(hits_folder)

####I commented this function bc it was too slow...
# def count_residues_in_chain(pdb_file, chain_id='A'):
#     # Create a PDB parser
#     parser = PDB.PDBParser(QUIET=True)

#     # Load the structure
#     structure = parser.get_structure('pdb', pdb_file)

#     # Initialize a counter for residues
#     residue_count = 0

#     # Iterate through each model in the structure
#     for model in structure:
#         # Iterate through each chain in the model
#         for chain in model:
#             # Check if the chain ID matches the desired chain
#             if chain.id == chain_id:
#                 # Count the number of residues in the chain
#                 residue_count += sum(1 for _ in chain.get_residues())

#     return residue_count

####This one is a bit messier and dirtier but does the job
def count_residues_in_chain(pdb_file, chain_id='A'):
    residue_count = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line.split()[4] == chain_id:
                residue_count += 1
    return residue_count/4


# Trim descriptions of silents to get original pdb designs names
def trim_substituted(text):
    file = text.split('_substituted')[0]+'.pdb'
    whole_path = 'output/'+file.split('_design')[0]+'/'+file
    return whole_path

# Function to merge CSV files in a directory (including subdirectories)
def merge_csv_files(directory):
    df_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.sc'):
                file_path = os.path.join(root, file)
                df = pd.read_table(file_path, sep=r'\s+', encoding='utf-8')
                df_list.append(df)
    if df_list:
        merged_df = pd.concat(df_list)
        merged_df['original_design'] = merged_df['description'].apply(trim_substituted)
        merged_df['length'] = merged_df['original_design'].apply(count_residues_in_chain)
        return merged_df
    else:
        return pd.DataFrame(columns=['plddt_binder','pae_interaction'])

# Function to track job status
def track_job_status(directory):
    job_status = []
    for root, dirs, files in os.walk(directory):
        for dir_name in dirs:
            job_path = os.path.join(root, dir_name)
            trj_path = os.path.join(job_path, "trj")  # Path to the "trj" directory

            job_record = {"job": dir_name, "status": "Waiting"}  # Default status

            # Check for RFD status
            if os.path.exists(trj_path) and os.path.isdir(trj_path):
                job_record["status"] = "RFD" 

            # Check for pMPNN status
            input_files = [f for f in os.listdir(job_path) if f.endswith('_input.silent.idx')]
            if input_files:
                job_record["status"] = "pMPNN"

            # Check for AF2 status
            af2_files = [f for f in os.listdir(job_path) if f.endswith('_input_out.silent.idx')]
            if af2_files:
                job_record["status"] = "AF2"

            # Check if the job has finished
            done_file = [f for f in os.listdir(job_path) if f == dir_name + '_done']
            sc_files = [f for f in os.listdir(job_path) if f.endswith('.sc')]

            if done_file:
                if sc_files:  # If there are *.sc files, consider it finished
                    job_record["status"] = "Finished"
                else:  # If there are no *.sc files, mark as FAILED
                    job_record["status"] = "FAILED"

            job_status.append(job_record)
        break  # Exit after the first iteration to limit depth
    
    for record in job_status:
        job_name = record["job"]
        numeric_part = int(re.search(r'\d+', job_name).group())  # Extract numeric part
        record["sort_key"] = numeric_part  # Add as sort_key

    # Sort by this new key before passing to DataTable
    status_df = pd.DataFrame(job_status).sort_values(by="sort_key")
    return status_df

# Function to update scatter plot, histograms, and row count
def update_scatter_plot(directory):
    merged_df = merge_csv_files(directory)
    if not merged_df.empty:
        filtered_df = merged_df[(merged_df['plddt_binder'] >= plddt_thres) & (merged_df['pae_interaction'] <= pae_thres)]
        row_count_text = f"Finished models: {len(merged_df)} | Number of hits: {len(filtered_df)} | Hit efficiency: {round(len(filtered_df)/len(merged_df)*100, 2)}%"
        
        scatter_plot = px.scatter(merged_df,
                                  y='plddt_binder',
                                  x='pae_interaction',
                                  color='length',
                                  marginal_x = "violin",
                                  marginal_y = "violin",
                                  render_mode = 'webgl',
                                  hover_data=['original_design','plddt_binder','pae_interaction','length']
                                  )
        scatter_plot.update_xaxes(range=[0,30], autorange="reversed")
        scatter_plot.update_yaxes(range=[40,100])
        scatter_plot.update_traces(marker=dict(opacity=0.7))
        # scatter_plot.add_hline(y=plddt_thres, line_dash="dot", line_color="grey")
        # scatter_plot.add_vline(x=pae_thres, line_dash="dot", line_color="grey")
        scatter_plot.update_layout(
            margin=dict(l=80, r=80, t=50, b=100),
            showlegend=True,
            legend=dict(x=1, y=1),
            shapes=[
                dict(type="rect", x0=0, y0=plddt_thres, x1=pae_thres, y1=100, line=dict(color="black"), opacity=0.5),
                # dict(type="line", xref="x", yref="paper", x0=10, y0=0, x1=10, y1=1, line=dict(color="grey", width=2, dash="dot")),
                # dict(type="line", xref="paper", yref="y", x0=0, y0=80, x1=1, y1=80, line=dict(color="grey", width=2, dash="dot")),
                ]
            )
        return scatter_plot, row_count_text
    else:
        return None, "No data available."

# Function to get best hit name
def get_hit_names(df, pae_interaction_threshold=pae_thres):
    filtered_df = df[df['pae_interaction'] <= pae_interaction_threshold]
    if not filtered_df.empty:
        sorted_df = filtered_df.sort_values(by='pae_interaction', ascending=True)
        hit_names = sorted_df['description'].tolist()
        return hit_names
    else:
        return ["No hits found under the specified conditions"]

# Function to get best hit path
def get_design_file_path_and_name(description, directory):
    match = re.match(r"(run_\d+)_design_(\d+_substituted).*", description)
    if not match:
        print("Invalid description format")
        # Returning None for both values if the format is invalid
        return None, None

    run_part, design_part = match.groups()
    # Constructing the file path using 'directory' which is now an absolute path
    data_path = os.path.join(directory, run_part)
    filename = f"{run_part}_design_{design_part}"

    return data_path, filename

# Styles
bt_style = {"align-items": "center", "background-color": "#F2F3F4", "border": "2px solid #000",
            "box-sizing": "border-box", "color": "#000", "cursor": "pointer", "display": "inline-flex",
            "font-family": "Helvetica", 'padding':'0.3em 1.2em', 'margin':'0.5em 0.3em 0.3em 0',
            "font-size": "0.9em", 'font-weight':'500', 'border-radius':'2em', 'text-align':'center',
            'transition':'all 0.2s'}
table_viewer_div_style={'display': 'flex', 'justifyContent': 'flex-start'}
table_div_style={'marginRight': '30px'}
table_style={'overflowX': 'auto','width': '100%', 'margin': '0 auto'}
table_cell_style={'minWidth': '150px', 'width': '150px', 'maxWidth': '150px', 'overflow': 'hidden', 'textOverflow': 'ellipsis', 'padding': '5px',}
title_style = {"margin-left": "15px", "margin-top": "15px", "margin-bottom": "0em", "color": "Black", "font-family" : "Helvetica", "font-size":"2.5em"}
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
app.layout = html.Div([
    html.H1("RFD Microruns watcher", style=title_style),
    html.Div(id='row-count'),
    html.Button('Make CSVs & Extract hits', id='execute-button', n_clicks=0, style=bt_style),
    dcc.Graph(id='scatter-plot', style={'width': '80vh', 'height': '70vh'}),
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
        html.Div([
            dcc.Dropdown(options='', value='', placeholder='Select a hit', id='ngl-molecule-dropdown'),
            dashbio.NglMoleculeViewer(id="ngl-molecule"),
            ]),
        ],style=table_viewer_div_style),
    dcc.Interval(
        id='interval-component',
        interval=180000,  # in milliseconds, update every 3 minutes
        n_intervals=0
    )
])

# PDB viewer
@callback(
    Output("ngl-molecule", 'data'),
    Output("ngl-molecule", "molStyles"),
    Input("ngl-molecule-dropdown", "value")
)

def return_molecule(value):
    if (value is None) or (value == ''):
        raise PreventUpdate
    else:
        data_path, filename = get_design_file_path_and_name(value, output_dir)
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
@app.callback(
    [Output('scatter-plot', 'figure'),
     Output('row-count', 'children'),
     Output('job-status-table', 'data'),
     Output('job-status-counts', 'children'),
     Output('ngl-molecule-dropdown', 'options')],
    [Input('interval-component', 'n_intervals'),
     Input('execute-button', 'n_clicks')])

def update_graph(n, n_clicks):
    directory='./output/'
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(directory)
    status_df = pd.DataFrame()
    status_df = track_job_status(directory)
    if not status_df.empty:
        status_df_records = status_df.to_dict('records')
    else:
        status_df_records = []

    ctx = dash.callback_context
    if ctx.triggered[0]['prop_id'] == 'execute-button.n_clicks':
        if not merged_df.empty:
            merged_df.iloc[:, 1:].to_csv('all_models.csv', index=False)
            filtered_df = merged_df[(merged_df['plddt_binder'] >= plddt_thres) & (merged_df['pae_interaction'] <= pae_thres)]
            if not filtered_df.iloc[:, 1:].empty:
                filtered_df.to_csv('hits.csv', index=False)
                for index, row in filtered_df.iterrows():
                    print('###############################\nEXTRACTING HIT\n###############################\n' + str(row))
                    run_number = re.search(r'run_\d+', row['description']).group()
                    command = f'silentextractspecific ../output/'+ run_number + '/' + run_number + '_input_out_af2.silent ' + row['description'] + ' > extraction.log'
                    subprocess.run(command, cwd=hits_folder, shell=True)
        scatter_plot, row_count_text = update_scatter_plot(directory)
    else:
        scatter_plot, row_count_text = update_scatter_plot(directory)

    # Job status counts
    status_counts = status_df['status'].value_counts().to_dict()
    status_counts_formatted = ", ".join([f"{status}: {count}" for status, count in status_counts.items()])
    job_status_counts_text = f"Job Status Counts: {status_counts_formatted}"

    # Dropdown HITS
    dropdown_options = get_hit_names(merged_df)

    return scatter_plot, row_count_text, status_df_records, job_status_counts_text, dropdown_options

# Run the app
if __name__ == '__main__':
    app.run_server(debug=debug_mode, dev_tools_hot_reload = False, use_reloader=True,
                   host=hostname, port=port_number)
