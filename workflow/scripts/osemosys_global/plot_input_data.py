import pandas as pd
import plotly as py
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import itertools
import os
import sys
import yaml
from osemosys_global.OPG_configuration import ConfigFile, ConfigPaths

pd.set_option('mode.chained_assignment', None)


def main():
    '''Creates figures of selected input data'''

    config_paths = ConfigPaths()
    config = ConfigFile('config')
    data_dir = config_paths.output_data_dir
    figure_dir = config_paths.output_dir

    # Check for output directory
    try:
        os.makedirs(os.path.join(figure_dir, 'figs'))
    except FileExistsError:
        pass

    for each_tech in ['SPV', 'WON', 'HYD']:
        plot_capfac(data_dir, each_tech, figure_dir)


def plot_capfac(data_dir, tech, figure_dir):
    df = pd.read_csv(os.path.join(data_dir,
                                  'CapacityFactor.csv'))

    #df = df.loc[df.TECHNOLOGY.str[6:9].isin(['IND',
    #                                         'BTN',
    #                                         'BGD',
    #                                         'NPL'])]
    df = df.loc[(df.TECHNOLOGY.str.startswith('PWR')) &
                (df.TECHNOLOGY.str[3:6].isin([tech]))]
    df.VALUE = df.VALUE.astype('float64')

    df.drop(['REGION'],
            axis=1,
            inplace=True)

    df = df.groupby(['TIMESLICE',
                     'TECHNOLOGY'],
                    as_index=False)['VALUE'].mean()

    fig = go.Figure(data=go.Heatmap(
        z=df['VALUE'],
        x=df['TIMESLICE'],
        y=df['TECHNOLOGY'],
        colorscale='inferno'
                        ),
                    )
    
    fig.update_layout(
        xaxis_nticks=len(df['TECHNOLOGY']),
        height=1400)
    # fig.update_xaxes(title_text='Hours')

    fig.write_image(os.path.join(figure_dir,
                                 'figs',
                                 'capfac_' +
                                 #each_node +
                                 ' ' +
                                 tech +
                                 '.jpg'))


if __name__ == '__main__':
    main()
