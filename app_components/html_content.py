from dash_core_components import Markdown, Upload
from dash_html_components import (
  Div, H3, P, A,
  Button as html_Button, Img
)
from dash_bootstrap_components import (
  Container,
  Row, Col,
  Navbar, NavbarBrand,
  Card, CardHeader, CardBody,
  Tabs, Tab,
  Button as dbc_Button
)
from dash import Dash


def create_app_skeleton(header: str, fasta_input_card: Card, info_card: Card, table_card: Card):
  return (
    Div(
      [
        header,
        Div(Markdown('''
                     If you use EffectorO in your research, please cite us:
                      
                     Nur, M. J., Wood, K. J. & Michelmore, R. W. EffectorO: motif-independent prediction of effectors in oomycete genomes using machine learning and lineage-specificity. Mol. Plant-Microbe Interact. (2023). doi:10.1094/MPMI-11-22-0236-TA

                     Also — if you are using EffectorO to predict oomycete effectors
                     we would love to hear from you! Please email Kelsey at
                     [klsywd@gmail.com](mailto:klsywd@gmail.com) or on Twitter
                     [@klsywd](https://twitter.com/klsywd) and let us know what
                     organism you are studying.
                     '''),
            style={'margin': '30px'}
          ),
        # Interval(id='interval-component',
        #              interval=2000,  # Refresh the memory usage every 2 seconds
        #              n_intervals=0),
        Container([Row([Col(fasta_input_card, md=6),
                        Col(info_card, md=6)]),
                   Row([Col(table_card)])],
			                fluid=True),
        Div(id="memory-usage"),
      ],
    )
  )


# title navbar

def create_title_navbar(app: Dash):
  genomecenter_icon = A(
      Img(src=app.get_asset_url("genomecenter_logo.png"), height='30px'),
      href="https://genomecenter.ucdavis.edu/"
    )
  
  github_button = dbc_Button(
    "View code on GitHub",
    outline=False,
    color="primary",
    href="https://github.com/mjnur/oomycete-effector-prediction",
    id="gh-link",
    style={"text-transform": "none"},
    )

  return (
    Navbar(
      Container(
        [
          A(genomecenter_icon, className="navbar-brand",
            href="https://michelmorelab.ucdavis.edu/effectoro"),
          NavbarBrand("EffectorO: Motif-Independent Oomycete Effector Prediction"),
          # A(
          #   [genomecenter_icon,
          #    " ", # space padding
          #    "EffectorO: Motif-Independent Oomycete Effector Prediction"],
          #   className="navbar-brand",
          #   href="https://michelmorelab.ucdavis.edu/effectoro",
          #   style={'text-decoration': 'none'}
          # ),
          github_button
        ],
        fluid=True
      ),
      style={'align': 'center'},
      color='dark',
      dark=True
    )
  )


# card creation functions

def create_fasta_input_card(max_byte_size):
  return (
    Card(
      [
        CardHeader(H3("Input FASTA protein file for EffectorO-ML")),
        CardBody(
          Row([
            Col([
              # upload button with max file size 10 MB
              Upload(
                html_Button(
                  "Select a FASTA file",
                  # TODO: edit style of input FASTA file button
                ),
                id="upload-data",
                max_size=max_byte_size,
                # TODO: edit style of upload
              )
            ]),
            Div([
              Markdown(
                f'''
                 **File Requirements:**

                 - FASTA-formatted file of predicted amino acid sequences
                 - File less than {max_byte_size / 10**6} MB
                 '''
              )
            ])
          ])
        )
        # TODO: work on this
        # inspirations from POOE:
        #   - implement ability to copy paste FASTA file input
        #   - specify job name?
      ],
      # TODO: style of FASTA input card
    )
  )

def create_info_card():
  return (
    Card(
      [
        CardHeader(
          Tabs(
            [Tab(label="Launch Instructions", tab_id="launch-tab"),
            Tab(label="Methodology", tab_id="method-tab")],
            id="info-card-tabs",
            card=True,
            active_tab="launch-tab",
            style={'cursor': 'pointer'}
          )
        ),
        CardBody(P(id="card-content", className="card-text"))
      ],
      # TODO: info card style
    )
  )

def create_table_card():
  return (
    Card(
      [
        CardHeader(H3("EffectorO-ML Prediction Table")),
        CardBody(Row(Col([Div(id="datatable")])))
      ]
    )
  )