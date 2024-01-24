import pandas as pd
from psutil import Process
from numpy import array, round
from base64 import b64decode
from Bio.SeqIO import parse
from io import StringIO

from dash.dependencies import Input, Output, State
from dash_table import DataTable
from dash_core_components import Markdown
from dash_html_components import Div

from app_components.get_average_features import get_average_features

def get_callbacks(app):
	def parse_contents(contents, trained_model):
		try:
			content_type, content_string = contents.split(',')
			decoded = b64decode(content_string)
			seqs_to_predict = parse(StringIO(decoded.decode('utf-8')),'fasta')
			prediction_map = {'0': "predicted non-effector", '1': "predicted effector"}

			seq_ids = []
			full_sequences = []
			for protein in seqs_to_predict:
				seq_ids.append(protein.id)
				full_sequences.append(protein.seq)
			seq_features = array([get_average_features(seq) for seq in full_sequences])
			# get predicted output
			predictions = trained_model.predict(seq_features)
			probabilities = trained_model.predict_proba(seq_features)
			meanings = array([prediction_map[pred] for pred in predictions])

			df = pd.DataFrame({"proteinID": seq_ids,
												"prediction": predictions,
												"probability": probabilities[:,1],
												"meaning": meanings
												})

			# round probabilities
			df['probability'] = round(df['probability'], 3)

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

		# TODO: use pd.DataFrame.rename(...) here
		table['Protein ID'] = table['proteinID']
		table['Probability'] = table['probability']
		table['Classification'] = table['prediction']
		table['Prediction'] = table['meaning']
		table = table[['Protein ID', 'Probability', 'Classification', 'Prediction']]

		formatted_table = DataTable(
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


		return Div([formatted_table])


	def get_memory_usage():
		return Process().memory_info().rss  # Get the resident set size (RSS) memory usage

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
			return Div([Markdown('''
													See [published paper]
													(https://apsjournals.apsnet.org/doi/10.1094/MPMI-11-22-0236-TA)
													for detailed information.
												
													This pipeline runs **EffectorO-ML**, a pre-trained
													machine-learning based Oomycete effector
													classifier, built from Random Forest models using
													biochemical amino acid characteristics as
													features.
													'''),])

		elif active_tab == "launch-tab":
			return Div([Markdown('''
													1. Click on the "**SELECT A FASTA FILE**"
													upload box.

													2. Upload a FASTA file of predicted amino
													acid sequences, for EffectorO-ML to analyze.

													3. View the sortable and filterable results
													in the datatable below.
													'''),])


	# callback to toggle the collapse on small screens
	@app.callback(
		Output("navbar-collapse", "is_open"),
		[Input("navbar-toggler", "n_clicks")],
		[State("navbar-collapse", "is_open")],
	)
	def toggle_navbar_collapse(n, is_open):
		return is_open if (n) else not is_open