import pandas as pd
pd.set_option('mode.chained_assignment', None)
import plotly.express as px
import matplotlib.pyplot as plt
import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from mpl_toolkits.basemap import Basemap
from typing import Dict
from osemosys_global.utils import read_csv
from osemosys_global.visualisation.utils import get_color_codes
from osemosys_global.visualisation.data import get_total_capacity_data, get_generation_annual_data, get_generation_ts_data
from osemosys_global.configuration import ConfigFile, ConfigPaths


def main():
    """Creates system level and country level graphs."""

    config_paths = ConfigPaths()
    config = ConfigFile('config')
    scenario_figs_dir = config_paths.scenario_figs_dir
    results_by_country = config.get('results_by_country')

    # Check for output directory 
    try:
        os.makedirs(scenario_figs_dir)
    except FileExistsError:
        pass
    
    input_data = read_csv(config_paths.scenario_data_dir)
    result_data = read_csv(config_paths.scenario_results_dir)
    
    # Get system level results 
    plot_total_capacity(result_data, scenario_figs_dir, country=None)
    plot_generation_annual(result_data, scenario_figs_dir, country=None)
    plot_generation_hourly(input_data, result_data, scenario_figs_dir, country=None)

    # If producing by country results, check for folder structure 
    if results_by_country:
        countries = config.get('geographic_scope')
        for country in countries:
            try:
                os.makedirs(os.path.join(scenario_figs_dir, country))
            except FileExistsError:
                pass
    
            plot_total_capacity(result_data, scenario_figs_dir, country=country)
            plot_generation_annual(result_data, scenario_figs_dir, country=country)
    
    # Creates transmission maps by year      
    years = [config.get('endYear')]
    for year in years:
        plot_transmission_capacity(year)
        plot_transmission_flow(year)

def plot_total_capacity(data: Dict[str,pd.DataFrame], save_dir: str, country:str = None):
    """Plots total capacity chart 
        
    Arguments:
        save_dir: str 
            Directory to save file 
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    """
    
    df = get_total_capacity_data(data, country=country)
    plot_colors = get_color_codes()

    if not country: # System level titles
        graph_title = 'Total System Capacity'
        legend_title = 'Powerplant'
    else: # Country level titles
        graph_title = f'{country} System Capacity'
        legend_title = 'Country-Powerplant'

    fig = px.bar(df,
                 x='YEAR',
                 y='VALUE',
                 color='LABEL',
                 color_discrete_map=plot_colors,
                 template='plotly_white',
                 labels={'YEAR': 'Year',
                         'VALUE': 'Gigawatts (GW)',
                         'LABEL': legend_title},
                 title=graph_title)
    
    fig.update_layout(
        font_family="Arial",
        font_size=14,
        legend_traceorder="reversed",
        title_x=0.5)
    fig['layout']['title']['font'] = dict(size=24)
    fig.update_traces(marker_line_width=0, opacity=0.8)

    if country:
        fig_file = os.path.join(save_dir, country, 'TotalCapacityAnnual.html')
    else:
        fig_file = os.path.join(save_dir, 'TotalCapacityAnnual.html')

    return fig.write_html(fig_file)

def plot_generation_annual(data: Dict[str,pd.DataFrame], save_dir: str, country:str = None):
    """Plots total annual generation
        
    Arguments:
        save_dir: str 
            Directory to save file 
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    """

    df = get_generation_annual_data(data, country=None)
    plot_colors = get_color_codes()

    if not country: # System level titles
        graph_title = 'Total System Generation'
        legend_title = 'Powerplant'
    else: # Country level titles
        graph_title = f'{country} System Generation'
        legend_title = 'Country-Powerplant'

    fig = px.bar(df,
                 x='YEAR',
                 y='VALUE',
                 color='LABEL',
                 color_discrete_map=plot_colors,
                 template='plotly_white',
                 labels={'YEAR': 'Year',
                         'VALUE': 'Petajoules (PJ)',
                         'LABEL': legend_title},
                 title=graph_title)
    fig.update_layout(
        font_family="Arial",
        font_size=14,
        legend_traceorder="reversed",
        title_x=0.5)
    fig['layout']['title']['font'] = dict(size=24)
    fig.update_traces(marker_line_width=0,
                      opacity=0.8)

    if country:
        fig_file = os.path.join(save_dir, country, 'GenerationAnnual.html')
    else:
        fig_file = os.path.join(save_dir, 'GenerationAnnual.html')

    return fig.write_html(fig_file)


def plot_generation_hourly(
    input_data: Dict[str,pd.DataFrame], 
    result_data: Dict[str,pd.DataFrame], 
    save_dir: str, 
    country:str = None
):
    """Plots total annual generation
        
    Arguments:
        save_dir: str 
            Directory to save file 
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    """

    df = get_generation_ts_data(input_data, result_data, country=country)
    plot_colors = get_color_codes()

    fig = px.area(df,
                  x='HOUR',
                  y=[x for
                     x in df.columns
                     if x not in ['MONTH', 'HOUR']
                    ],
                  title='System Hourly Generation',
                  facet_col='MONTH',
                  facet_col_spacing=0.005,
                  color_discrete_map=plot_colors,
                  animation_frame='YEAR',
                  template='seaborn+plotly_white',
                  labels={
                      "variable":"",})
    fig.update_layout(
        legend_traceorder="reversed",
        title_x=0.5,
        yaxis_title = 'Gigawatt-Hours (GWh)')
    fig['layout']['title']['font'] = dict(size=24)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    '''
    for axis in fig.layout:
        if type(fig.layout[axis]) == go.layout.XAxis:
            fig.layout[axis].title.text = ''

        fig['layout']['yaxis']['title']['text']='Gigawatts (GW)'
        fig['layout']['yaxis3']['title']['text']=''
        fig['layout']['xaxis']['title']['text']=''
        fig['layout']['xaxis7']['title']['text']='Hours'
    '''
    if country:
        fig_file = os.path.join(save_dir, country, 'GenerationHourly.html')
    else:
        fig_file = os.path.join(save_dir, 'GenerationHourly.html')
    return fig.write_html(fig_file)
    

def midpoint(x1, y1, x2, y2):
    return ((x1 + x2)/2, (y1 + y2)/2)

def plot_transmission_capacity(year):

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_figs_dir = config_paths.scenario_figs_dir
    scenario_results_dir = config_paths.scenario_results_dir
    input_data_dir = config_paths.input_data_dir
    
    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'TotalCapacityAnnual.csv'
                                  )
                     )
    
    df_centerpoints = pd.read_excel(os.path.join(input_data_dir,
                                                'Costs Line expansion.xlsx',
                                                ), sheet_name = 'Centerpoints'
                                    )
    
    df_centerpoints['Node'] = df_centerpoints['Node'].str.split('-', n = 1).str[1]
    df_centerpoints['Node'] = np.where(df_centerpoints['Node'].str.len()==3, 
                                       df_centerpoints['Node'] + 'XX', 
                                       df_centerpoints['Node'].str.replace('-', '')
                                       )
    
    df_centerpoints.set_index('Node', inplace = True)
    
    df = df.loc[df['TECHNOLOGY'].str.startswith('TRN')]
    df.VALUE = df.VALUE.astype('float64')
    df['FROM'], df['TO'] = df.TECHNOLOGY.str[3:8], df.TECHNOLOGY.str[8:]
    
    df = df.merge(df_centerpoints[['Latitude', 'Longitude']], left_on = 'FROM', right_index = True)
    df = df.merge(df_centerpoints[['Latitude', 'Longitude']], left_on = 'TO', right_index = True, suffixes = ('_FROM', '_TO'))
    df = df.groupby(['TECHNOLOGY', 'YEAR', 'FROM', 'TO', 'Longitude_FROM', 'Latitude_FROM', 'Longitude_TO', 'Latitude_TO'],
                    as_index=False)['VALUE'].sum()
    
    scaler = MinMaxScaler()
    maxlinewidth = 3
    df['line_width'] = (scaler.fit_transform(df[['VALUE']]) * maxlinewidth).round(1)
    df['line_width'].where(df['line_width'] >= 0.3, 0.3, inplace = True)
    
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    
    dfyear = df.loc[df['YEAR'] == year]
    
    maplatmin = int(df[['Latitude_FROM', 'Latitude_TO']].min().min() +-6)
    maplatmax = int(df[['Latitude_FROM', 'Latitude_TO']].max().max() +6)
    maplonmin = int(df[['Longitude_FROM', 'Longitude_TO']].min().min() +-6)
    maplonmax = int(df[['Longitude_FROM', 'Longitude_TO']].max().max() +6)
    
    # create new figure
    fig, ax = plt.subplots()
    #ax.set_title(f'Total Transmission Capacity in {year} (GW)', fontsize = '6')
    npa = np.array #basemap has limited input data options, doesn't read floats
    
    # setup mercator map projection.
    m = Basemap(projection='cyl',llcrnrlat=maplatmin,urcrnrlat=maplatmax,\
                llcrnrlon=maplonmin,urcrnrlon=maplonmax,resolution='c'
                )
    
    # generates all lines and map text
    for y in dfyear.index.unique():
        Map_Filter = dfyear.loc[(dfyear.index == y)]
        x1,y1 = npa(Map_Filter['Longitude_FROM']) , npa(Map_Filter['Latitude_FROM'])
        x2,y2 = npa(Map_Filter['Longitude_TO']) , npa(Map_Filter['Latitude_TO'])
        
        midco = midpoint(x1, y1, x2, y2)
        
        m.drawgreatcircle(x1 , y1 , x2, y2, linewidth = npa(Map_Filter['line_width']) , color = 'crimson')       
        m.drawcoastlines(color = 'black', linewidth = 0.3)
        m.drawmapboundary(fill_color='#46bcec')
        m.fillcontinents(color = 'lightgrey')
        m.drawcountries(color="black", linewidth = 0.3)
       # plt.text(x1, y1, Map_Filter['FROM'].iloc[0], fontsize = 2, fontweight= 'bold', ha = 'center', va = 'top', alpha = .8)
       # plt.text(x2, y2, Map_Filter['TO'].iloc[0], fontsize = 2, fontweight= 'bold', ha = 'center', va = 'top', alpha = .8)
        plt.text(midco[0], midco[1], int(npa(Map_Filter['VALUE'].iloc[0])), fontsize = 2.5, fontweight= 'bold', ha = 'right', va = 'top', alpha = .8)

    return fig.savefig(os.path.join(scenario_figs_dir,
                         f'TransmissionCapacity{year}.jpg'
                         ), dpi = 500, bbox_inches = 'tight'
            )

def plot_transmission_flow(year):

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_figs_dir = config_paths.scenario_figs_dir
    scenario_results_dir = config_paths.scenario_results_dir
    input_data_dir = config_paths.input_data_dir
    
    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'ProductionByTechnologyAnnual.csv'
                                  )
                     )
    
    df_centerpoints = pd.read_excel(os.path.join(input_data_dir,
                                                'Costs Line expansion.xlsx',
                                                ), sheet_name = 'Centerpoints'
                                    )
    
    df_centerpoints['Node'] = df_centerpoints['Node'].str.split('-', n = 1).str[1]
    df_centerpoints['Node'] = np.where(df_centerpoints['Node'].str.len()==3, 
                                       df_centerpoints['Node'] + 'XX', 
                                       df_centerpoints['Node'].str.replace('-', '')
                                       )
    
    df_centerpoints.set_index('Node', inplace = True)
    
    df = df.loc[df['TECHNOLOGY'].str.startswith('TRN')]
    df.VALUE = df.VALUE.astype('float64')
    df['FROM'], df['TO'] = df.TECHNOLOGY.str[3:8], df.TECHNOLOGY.str[8:]
    
    df = df.merge(df_centerpoints[['Latitude', 'Longitude']], left_on = 'FROM', right_index = True)
    df = df.merge(df_centerpoints[['Latitude', 'Longitude']], left_on = 'TO', right_index = True, suffixes = ('_FROM', '_TO'))
    df = df.groupby(['TECHNOLOGY', 'YEAR', 'FROM', 'TO', 'Longitude_FROM', 'Latitude_FROM', 'Longitude_TO', 'Latitude_TO'],
                    as_index=False)['VALUE'].sum()
    
    scaler = MinMaxScaler()
    maxlinewidth = 3
    df['line_width'] = (scaler.fit_transform(df[['VALUE']]) * maxlinewidth).round(1)
    df['line_width'].where(df['line_width'] >= 0.3, 0.3, inplace = True)
    
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    
    dfyear = df.loc[df['YEAR'] == year]
    
    maplatmin = int(df[['Latitude_FROM', 'Latitude_TO']].min().min() +-6)
    maplatmax = int(df[['Latitude_FROM', 'Latitude_TO']].max().max() +6)
    maplonmin = int(df[['Longitude_FROM', 'Longitude_TO']].min().min() +-6)
    maplonmax = int(df[['Longitude_FROM', 'Longitude_TO']].max().max() +6)
    
    # create new figure
    fig, ax = plt.subplots()
    #ax.set_title(f'Total Transmission Flow in {year} (PJ)', fontsize = '6')
    npa = np.array #basemap has limited input data options, doesn't read floats
    
    # setup mercator map projection.
    m = Basemap(projection='cyl',llcrnrlat=maplatmin,urcrnrlat=maplatmax,\
                llcrnrlon=maplonmin,urcrnrlon=maplonmax,resolution='c'
                )
    
    # generates all lines and map text
    for y in dfyear.index.unique():
        Map_Filter = dfyear.loc[(dfyear.index == y)]
        x1,y1 = npa(Map_Filter['Longitude_FROM']) , npa(Map_Filter['Latitude_FROM'])
        x2,y2 = npa(Map_Filter['Longitude_TO']) , npa(Map_Filter['Latitude_TO'])
        
        midco = midpoint(x1, y1, x2, y2)
        
        m.drawgreatcircle(x1 , y1 , x2, y2, linewidth = npa(Map_Filter['line_width']) , color = 'crimson')       
        m.drawcoastlines(color = 'black', linewidth = 0.3)
        m.drawmapboundary(fill_color='#46bcec')
        m.fillcontinents(color = 'lightgrey')
        m.drawcountries(color="black", linewidth = 0.3)
       # plt.text(x1, y1, Map_Filter['FROM'].iloc[0], fontsize = 2, fontweight= 'bold', ha = 'center', va = 'top', alpha = .8)
       # plt.text(x2, y2, Map_Filter['TO'].iloc[0], fontsize = 2, fontweight= 'bold', ha = 'center', va = 'top', alpha = .8)
        plt.text(midco[0], midco[1], int(npa(Map_Filter['VALUE'].iloc[0])), fontsize = 2.5, fontweight= 'bold', ha = 'right', va = 'top', alpha = .8)

    return fig.savefig(os.path.join(scenario_figs_dir,
                         f'TransmissionFlow{year}.jpg'
                         ), dpi = 500, bbox_inches = 'tight'
            )

if __name__ == '__main__':
    main()
