import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

def create_app_skeleton(header, fasta_input_card, info_card, table_card):
  return (
    html.Div(
      [
        header,
        dcc.Markdown('''If you use EffectorO in your research, please cite us:'''),
        html.Blockquote('Nur, M. J., Wood, K. J. & Michelmore, R. W. EffectorO: motif-independent prediction of effectors in oomycete genomes using machine learning and lineage-specificity. Mol. Plant-Microbe Interact. (2023). doi:10.1094/MPMI-11-22-0236-TA'),
        dcc.Markdown('''
                     Also — if you are using EffectorO to predict oomycete effectors
                     we would love to hear from you! Please email Kelsey at
                     [klsywd@gmail.com](mailto:klsywd@gmail.com) or on Twitter
                     [@klsywd](https://twitter.com/klsywd) and let us know what
                     organism you are studying.
                     '''),
        dcc.Interval(id='interval-component',
                     interval=2000,  # Refresh the memory usage every 2 seconds
                     n_intervals=0),
        dbc.Container([dbc.Row([dbc.Col(fasta_input_card, md=6),
								                dbc.Col(info_card, md=6)]),
			                 dbc.Row([dbc.Col(table_card)])],
			                fluid=True),
        html.Div(id="memory-usage"),
      ]

      # TODO: add margin fix + style
    )
  )