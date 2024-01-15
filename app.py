from dash import Dash
import dash_table
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

import numpy as np
import pandas as pd
import pickle
from Bio import SeqIO
import io
from scripts.machine_learning_classification.scripts import FEAT as FEAT


from psutil import Process
import base64

import app_components.html_content as html_content
import app_components.callback_functions as callback_functions

# Set up the app
external_stylesheets = [dbc.themes.BOOTSTRAP, "assets/object_properties_style.css"]
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)
server = app.server

# set up ML model
trained_model = pickle.load(open('scripts/machine_learning_classification/trained_models/RF_88_best.sav', 'rb'))


def get_memory_usage():
	return Process().memory_info().rss  # Get the resident set size (RSS) memory usage

def get_average_features(sequence, df=0, protID=0):
	'''
	Method: Fills feature columns in a dataframe with net averages

	Input:

		- sequence: amino acid string
		- df: dataframe
		- protID: sequence identifier
	'''

	# get net cumulative sums
	#print("calculating for", sequence)
	length = min(len(sequence), 900)

	gravy = hydrophobicity = exposed = disorder = bulkiness = interface = 0.0
	for ind, aa in enumerate(sequence):
		if ind == length: break

		if aa.upper() in FEAT.INTERFACE_DIC:
			gravy += FEAT.GRAVY_DIC[aa.upper()]
			hydrophobicity += FEAT.HYDRO_DIC[aa.upper()]
			exposed += FEAT.EXPOSED_DIC[aa.upper()]
			disorder += FEAT.DISORDER_DIC[aa.upper()]
			bulkiness += FEAT.BULKY_DIC[aa.upper()]
			interface += FEAT.INTERFACE_DIC[aa.upper()]

	# return averages
	return [gravy						/ length,
				  hydrophobicity 	/ length,
					exposed 				/ length,
					disorder 				/ length,
					bulkiness 			/ length,
					interface 			/ length]


# file upload function
def parse_contents(contents, filename):
	try:
		content_type, content_string = contents.split(',')
		decoded = base64.b64decode(content_string)
		seqs_to_predict = SeqIO.parse(io.StringIO(decoded.decode('utf-8')),'fasta')
		prediction_map = {'0': "predicted non-effector", '1': "predicted effector"}

		seq_ids = []
		full_sequences = []
		for protein in seqs_to_predict:
			seq_ids.append(protein.id)
			full_sequences.append(protein.seq)
		seq_features = np.array([get_average_features(seq) for seq in full_sequences])
		# get predicted output
		predictions = trained_model.predict(seq_features)
		probabilities = trained_model.predict_proba(seq_features)
		meanings = np.array([prediction_map[pred] for pred in predictions])

		df = pd.DataFrame({"proteinID": seq_ids,
											 "prediction": predictions,
											 "probability": probabilities[:,1],
											 "meaning": meanings
											})

		# round probabilities
		df['probability'] = np.round(df['probability'], 3)

	except Exception as e:
		print(e)
		df = pd.DataFrame({"proteinID": [],
											 "prediction": [],
											 "probability": [],
											 "meaning": []
											 })
	return df

@app.callback(Output('datatable', 'children'),
              [Input('upload-data', 'contents'),
               Input('upload-data', 'filename')])
def get_new_datatable(contents, filename):
	table = parse_contents(contents, filename)

	table['Protein ID'] = table['proteinID']
	table['Probability'] = table['probability']
	table['Classification'] = table['prediction']
	table['Prediction'] = table['meaning']
	table = table[['Protein ID', 'Probability', 'Classification', 'Prediction']]

	formatted_table = dash_table.DataTable(
		columns=[{"name": i, "id": i} for i in table.columns],
		data=table.to_dict("records"),
		tooltip_header={
			'Protein ID': 'Protein ID obtained from fasta file',
			'Classification': 'Binary classification value of non-effector (0) and effector(1), \
												 using a cutoff of 0.5',
			'Probability': 'Random Forest predicted probability of being an effector',
			'Prediction': 'Classification value in a readable format'
		},
		style_header={
			"textDecoration": "underline",
			"textDecorationStyle": "dotted",
			'backgroundColor': 'white',
			'fontWeight': 'bold',
			"font-size": "16px"
		},
		tooltip_delay=0,
		page_size=18,
		export_format="csv",
		sort_action='native',
		filter_action='native',
		tooltip_duration=None,
		style_table={"overflowY": "scroll"},
		fixed_rows={"headers": False, "data": 0},
		style_cell={"width": "85px",
								"font-size": "16px",
								"fontFamily": 'Sans-Serif'
		},
	)


	return html.Div([formatted_table])
	return parse_contents(contents, filename).to_html()

# Buttons
# button_howto = dbc.Button(
# 	"View code on GitHub",
# 	outline=False,
# 	color="primary",
# 	href="https://github.com/mjnur/oomycete-effector-prediction",
# 	id="gh-link",
# 	style={"text-transform": "none"},
# )

# Define Header Layout
header = html_content.create_title_navbar(app)
# header = dbc.Navbar(
# 	dbc.Container([
# 		dbc.Row(
# 			[
# 				dbc.Col(html.A(
# 						html.Img(
# 							src=app.get_asset_url("genomecenter_logo.png"),
# 							height="30px",
# 						),
# 						href="https://bremia.ucdavis.edu",
# 				)),
# 				dbc.Col(dbc.NavbarBrand("EffectorO: Motif-Independent Oomycete Effector Prediction")),
# 			],
# 			align="center",
# 		),

# 		dbc.Row(
# 			dbc.Col([
# 					dbc.NavbarToggler(id="navbar-toggler"),
# 					dbc.Collapse(
# 						dbc.Nav(
# 							[dbc.NavItem(button_howto)],
# 							className="ml-auto",
# 							navbar=True,
# 						),
# 						id="navbar-collapse",
# 						navbar=True,
# 					)
# 			]),
# 			align="center",
# 		),
# 	],
# 	fluid=True,
# 	),
# 	color="dark",
# 	dark=True,
# )

# Define Cards
fasta_input_card = html_content.create_fasta_input_card(10000000)
# fasta_input_card = dbc.Card(
# 	[
# 		dbc.CardHeader(html.H3("Input FASTA protein file for EffectorO-ML")),
# 		dbc.CardBody(
# 			dbc.Row([
# 				dbc.Col([
# 					# create upload button, with max file size of 10MB
# 					dcc.Upload(html.Button('Select a fasta file',
# 																	style={'width': "100%",
# 																				 'borderStyle': 'dashed',
# 																				 'font-family': 'sans-serif',
# 																				 'borderRadius': '5px',
# 																				 'borderWidth': '2px',}),
# 										 id="upload-data",
# 										 max_size=20000000,
# 										 style={"width": "100%"}),
# 				], width=12, style={"padding-left": 0}),
# 				html.Div([
# 					dcc.Markdown('''**File requirements**:''', style={'padding-top': '10px'}),
# 					dcc.Markdown('''
# 											 - FASTA-formatted file of predicted amino acid sequences
# 											 - File size less than 20MB
# 											 ''')
# 				]),
# 			]),
# 		),
# 	],
# 	style={'height': '100%'}
# )

tab_label_style={"font-size": "1.3rem"}

info_card = html_content.create_info_card()
# info_card = dbc.Card(
# 	[
# 		dbc.CardHeader(
# 			dbc.Tabs(
# 				[
# 					dbc.Tab(label="Launch Instructions", tab_id="launch-tab",
# 									label_style=tab_label_style),
# 					dbc.Tab(label="Methodology", tab_id="method-tab",
# 									label_style=tab_label_style),
# 				],
# 				id="info-card-tabs",
# 				card=True,
# 				active_tab="launch-tab",
# 			)
# 		),
# 		dbc.CardBody(html.P(id="card-content", className="card-text")),
# 	], style={'height': '100%'}
# )

# table_card = dbc.Card([
# 		dbc.CardHeader(html.H2("EffectorO-ML Prediction Table")),
# 		dbc.CardBody(dbc.Row(dbc.Col([html.Div(id='datatable')]))),
# ])
table_card = html_content.create_table_card()

@app.callback(
    Output('memory-usage', 'children'),
    [Input('interval-component', 'n_intervals')]
)
def update_memory_usage():
    memory_usage = get_memory_usage()
    return f"Memory Usage: {memory_usage / 1024 / 1024:.2f} MB"

@app.callback(
    Output("card-content", "children"),
		[Input("info-card-tabs", "active_tab")]
)
def tab_content(active_tab):
	if active_tab == "method-tab":
		return html.Div([dcc.Markdown('''
																	See [published paper]
																	(https://apsjournals.apsnet.org/doi/10.1094/MPMI-11-22-0236-TA)
																	for detailed information.
																	'''),
										 dcc.Markdown('''
																	This pipeline runs **EffectorO-ML**, a pre-trained
																	machine-learning based Oomycete effector
																	classifier, built from Random Forest models using
																	biochemical amino acid characteristics as
																	features.
																	'''),])

	elif active_tab == "launch-tab":
		return html.Div([dcc.Markdown('''
																	1. Click on the "**SELECT A FASTA FILE**"
																	upload box.
																	'''),
										 dcc.Markdown('''
																	2. Upload a FASTA file of predicted amino
																	acid sequences, for EffectorO-ML to analyze.
																	'''),
										 dcc.Markdown('''
																	3. View the sortable and filterable results
																	in the datatable below.
																	'''),])


# create app
app.layout = html_content.create_app_skeleton(header, fasta_input_card, info_card, table_card)

# initialize callback functions

# we use a callback to toggle the collapse on small screens
@app.callback(
	Output("navbar-collapse", "is_open"),
	[Input("navbar-toggler", "n_clicks")],
	[State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
	return is_open if (n) else not is_open


if __name__ == "__main__":
	app.run_server(debug=True)
