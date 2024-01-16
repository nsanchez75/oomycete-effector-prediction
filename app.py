import pickle

from dash import Dash
import dash_bootstrap_components as dbc

import app_components.html_content as html_content

from app_components.callback_functions import get_callbacks

# set up the app
external_stylesheets = [dbc.themes.BOOTSTRAP, "assets/object_properties_style.css"]
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)
server = app.server

# set up ML model
trained_model = pickle.load(open("machine_learning_classification/trained_models/RF_88_best.sav", 'rb'))

# define HTML contents
header = html_content.create_title_navbar(app)
fasta_input_card = html_content.create_fasta_input_card(10000000)
info_card = html_content.create_info_card()
table_card = html_content.create_table_card()

# create app
app.layout = html_content.create_app_skeleton(header, fasta_input_card, info_card, table_card)

# import callback functions after app had been initialized
get_callbacks(app)


if __name__ == "__main__":
	app.run_server(debug=True)
