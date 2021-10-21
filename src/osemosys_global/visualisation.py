import pandas as pd
import plotly as py
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import os
import sys

input_folder = '../../osemosys_global_model/OsemosysGlobal/results'
output_folder = '../../osemosys_global_model/OsemosysGlobal/figures'

try:
    os.makedirs(output_folder)
except FileExistsError:
    pass


def plot_totalcapacity():
    df = pd.read_csv(os.path.join(input_folder,
                                  'TotalCapacityAnnual.csv'
                                  )
                     )
    df.VALUE = df.VALUE.astype('float64')
    # df.TECHNOLOGY = df.TECHNOLOGY.replace(tech_names)
    df['TYPE'] = df.TECHNOLOGY.str[3:6]
    df.TYPE = df.TYPE[df.TYPE != 'TRN']
    df['COUNTRY'] = df.TECHNOLOGY.str[6:9]
    df['LABEL'] = df['COUNTRY'] + '-' + df['TYPE']
    df = df.groupby(['LABEL','YEAR'], as_index=False)['VALUE'].sum()

    fig = px.bar(df,
                 x='YEAR',
                 y='VALUE',
                 color='LABEL',
                 # color_discrete_map=color_dict,
                 template='plotly_white',
                 labels={'YEAR':'Year',
                     'VALUE':'Gigawatts (GW)',
                     'LABEL':'Country-Powerplant'})
    # fig.update_xaxes(type='category')
    fig.update_layout(
        font_family="Arial",
        font_size=14)
    fig.update_traces(marker_line_width=0, opacity=0.8)

    return fig.write_html(os.path.join(output_folder,
                                       'TotalCapacityAnnual.html'
                                       )
                          )


def plot_generationannual():
    df = pd.read_csv(os.path.join(input_folder,
                                  'ProductionByTechnologyAnnual.csv'
                                  )
                     )
    df.VALUE = df.VALUE.astype('float64')
    df['TYPE'] = df.TECHNOLOGY.str[3:6]
    df.TYPE = df.TYPE[df.TYPE != 'TRN']
    df['COUNTRY'] = df.TECHNOLOGY.str[6:9]
    df['LABEL'] = df['COUNTRY'] + '-' + df['TYPE']
    df = df.groupby(['LABEL','YEAR'], as_index=False)['VALUE'].sum()

    fig = px.bar(df,
                 x='YEAR',
                 y='VALUE',
                 color='LABEL',
                 # color_discrete_map=color_dict,
                 template='plotly_white',
                 labels={'YEAR':'Year',
                     'VALUE':'Gigawatts (GW)',
                     'LABEL':'Country-Powerplant'})
    # fig.update_xaxes(type='category')
    fig.update_layout(
        font_family="Arial",
        font_size=14)
    fig.update_traces(marker_line_width=0, opacity=0.8)

    return fig.write_html(os.path.join(output_folder,
                                        'GenerationAnnual.html'
                                        )
                           )
plot_totalcapacity()
plot_generationannual()

