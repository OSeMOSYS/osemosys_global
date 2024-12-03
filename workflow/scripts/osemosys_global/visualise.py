import pandas as pd
pd.set_option('mode.chained_assignment', None)
import plotly.express as px
import os
import sys 
from typing import List
from sklearn.preprocessing import MinMaxScaler
from typing import Dict
from osemosys_global.utils import read_csv, filter_transmission_techs
from osemosys_global.visualisation.utils import(get_color_codes, get_map, plot_map_trn_line, plot_map_text)
from osemosys_global.visualisation.data import(get_total_capacity_data, get_generation_annual_data, 
                                               get_generation_ts_data)
import osemosys_global.constants as constants
from configuration import ConfigFile, ConfigPaths


def main(
    input_data: pd.DataFrame,
    result_data: pd.DataFrame,
    centerpoints: pd.DataFrame,    
    custom_nodes: List[str],
    custom_nodes_centerpoints : pd.DataFrame,
    scenario_figs_dir: str,
    countries: List[str],
    results_by_country: bool = True,
    years: List[int] = [2050],
):
    """Creates system level and country level graphs."""

    # Check for output directory 
    try:
        os.makedirs(scenario_figs_dir)
    except FileExistsError:
        pass

    # Get system level results 
    plot_total_capacity(result_data, scenario_figs_dir, country=None)
    plot_generation_annual(result_data, scenario_figs_dir, country=None)
    plot_generation_hourly(input_data, result_data, scenario_figs_dir, country=None)

    # If producing by country results, check for folder structure 
    if results_by_country:
        for country in countries:
            try:
                os.makedirs(os.path.join(scenario_figs_dir, country))
            except FileExistsError:
                pass
    
            plot_total_capacity(result_data, scenario_figs_dir, country=country)
            plot_generation_annual(result_data, scenario_figs_dir, country=country)
    
    # Creates transmission maps by year      
    for year in years:
        plot_transmission_capacity(custom_nodes, centerpoints, custom_nodes_centerpoints, 
                                   result_data, scenario_figs_dir, year)
        
        plot_transmission_flow(custom_nodes, centerpoints, custom_nodes_centerpoints, 
                                   result_data, scenario_figs_dir, year)

def plot_total_capacity(data: Dict[str,pd.DataFrame], save_dir: str, country:str = None) -> None:
    """Plots total capacity chart 
        
    Arguments:
        data: Dict[str,pd.DataFrame]
            Result data
        save_dir: str 
            Directory to save file 
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    """

    df = get_total_capacity_data(data, country=country)
    plot_colors = constants.COLORS

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

def plot_generation_annual(data: Dict[str,pd.DataFrame], save_dir: str, country:str = None) -> None:
    """Plots total annual generation
        
    Arguments:
        save_dir: str 
            Directory to save file 
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    """

    df = get_generation_annual_data(data, country=country)
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
) -> None:
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

    if country:
        fig_file = os.path.join(save_dir, country, 'GenerationHourly.html')
    else:
        fig_file = os.path.join(save_dir, 'GenerationHourly.html')
    return fig.write_html(fig_file)
    
def midpoint(x1, y1, x2, y2):
    return ((x1 + x2)/2, (y1 + y2)/2)

def plot_transmission_capacity(
    custom_nodes: list,
    df_centerpoints: pd.DataFrame(),
    df_custom_nodes_centerpoints: pd.DataFrame(),
    result_data: Dict[str,pd.DataFrame], 
    save_dir: str, 
    year:int
) -> None:
    
    # get result data
    total_cap_annual = result_data["TotalCapacityAnnual"]
    trn = filter_transmission_techs(total_cap_annual)
    
    if trn.empty:
        return
    
    trn.VALUE = trn.VALUE.astype('float64')
    trn['FROM'], trn['TO'] = trn.TECHNOLOGY.str[3:8], trn.TECHNOLOGY.str[8:]
    
    # get node data 
    df_centerpoints = df_centerpoints.rename(columns = {'lat' : 'LATITUDE', 'long' : 'LONGITUDE'})
    df_centerpoints.set_index('region', inplace = True)
    
    if custom_nodes:            
        df_custom_nodes_centerpoints = df_custom_nodes_centerpoints.rename(columns = {'lat' : 'LATITUDE', 
                                                                                      'long' : 'LONGITUDE'})
        df_custom_nodes_centerpoints.set_index('region', inplace = True)
        df_centerpoints = pd.concat([df_centerpoints, df_custom_nodes_centerpoints])
        
        custom_centerpoints_missing = [val for val in custom_nodes if val not in list(df_centerpoints.index)]

        if custom_centerpoints_missing:
            print(f'{custom_centerpoints_missing} defined as custom nodes but not added to custom node centerpoints. Custom node transmission lines will not be plotted!')
                     
    # merge result data with node data 
    trn = trn.merge(df_centerpoints[['LATITUDE', 'LONGITUDE']], left_on = 'FROM', right_index = True)
    trn = trn.merge(df_centerpoints[['LATITUDE', 'LONGITUDE']], left_on = 'TO', right_index = True, suffixes = ('_FROM', '_TO'))
    trn = trn.groupby(['TECHNOLOGY', 'YEAR', 'FROM', 'TO', 'LONGITUDE_FROM', 'LATITUDE_FROM', 'LONGITUDE_TO', 'LATITUDE_TO'],
                    as_index=False)['VALUE'].sum()
    
    if trn.empty:
        return
    
    # assign line widths based on result data 
    scaler = MinMaxScaler()
    maxlinewidth = 3
    trn['line_width'] = (scaler.fit_transform(trn[['VALUE']]) * maxlinewidth).round(1)
    trn['line_width'].where(trn['line_width'] >= 0.3, 0.3, inplace = True)
    
    # get unique nodes to plot 
    nodes_to_plot = {}
    unique_nodes = set(trn["FROM"].to_list() + trn["TO"].to_list())
    for node in unique_nodes:
        nodes_to_plot[node] = [
            df_centerpoints.loc[node]["LATITUDE"],
            df_centerpoints.loc[node]["LONGITUDE"],
        ]
    
    # get plotting extent 
    lat_min = int(trn[['LATITUDE_FROM', 'LATITUDE_TO']].min().min() +-6)
    lat_max = int(trn[['LATITUDE_FROM', 'LATITUDE_TO']].max().max() +6)
    lon_min = int(trn[['LONGITUDE_FROM', 'LONGITUDE_TO']].min().min() +-6)
    lon_max = int(trn[['LONGITUDE_FROM', 'LONGITUDE_TO']].max().max() +6)
    extent = [lon_min, lon_max, lat_min, lat_max]
    
    # create new figure
    fig, ax = get_map(extent=extent)
    
    # generates all lines and map text
    df_year = trn.loc[trn['YEAR'] == int(year)]
    for y in df_year.index.unique():
        
        # get data to plot
        map_filter = df_year.loc[(df_year.index == y)]
        x1,y1 = map_filter['LONGITUDE_FROM'].iloc[0] , map_filter['LATITUDE_FROM'].iloc[0]
        x2,y2 = map_filter['LONGITUDE_TO'].iloc[0] , map_filter['LATITUDE_TO'].iloc[0]
        midco = midpoint(x1, y1, x2, y2)
        
        # add transmission line
        plot_map_trn_line(
            ax=ax,
            x1=x1,
            y1=y1,
            x2=x2,
            y2=y2,
            width=map_filter['line_width'].iloc[0]
        )
        
        # add capacity label
        plot_map_text(
            ax=ax,
            x=midco[0],
            y=midco[1],
            text=int(map_filter['VALUE'].iloc[0])
        )
        
    # add node labels 
    for _, node in enumerate(nodes_to_plot):
        plot_map_text(
            ax=ax,
            x=nodes_to_plot[node][1],
            y=nodes_to_plot[node][0],
            text=node,
        )
        
    ax.set_title(f'Transmission Capacity in {year} (GW)', fontsize = '6')
    return fig.savefig(
        os.path.join(save_dir, f'TransmissionCapacity{year}.jpg'), 
        dpi = 500, 
        bbox_inches="tight"
    )

def plot_transmission_flow(
    custom_nodes: list,
    df_centerpoints: pd.DataFrame(),
    df_custom_nodes_centerpoints: pd.DataFrame(),
    result_data: Dict[str,pd.DataFrame], 
    save_dir: str, 
    year:int
) -> None:

    # get result data
    production_annual = result_data["ProductionByTechnologyAnnual"]
    prd = filter_transmission_techs(production_annual)
    
    if prd.empty:
        return
    
    prd.VALUE = prd.VALUE.astype('float64')
    prd['FROM'], prd['TO'] = prd.TECHNOLOGY.str[3:8], prd.TECHNOLOGY.str[8:]
    
    # get node data 
    df_centerpoints = df_centerpoints.rename(columns = {'lat' : 'LATITUDE', 'long' : 'LONGITUDE'})
    df_centerpoints.set_index('region', inplace = True)
    
    if custom_nodes:
        df_custom_nodes_centerpoints = df_custom_nodes_centerpoints.rename(columns = {'lat' : 'LATITUDE', 
                                                                                      'long' : 'LONGITUDE'})
        df_custom_nodes_centerpoints.set_index('region', inplace = True)
        df_centerpoints = pd.concat([df_centerpoints, df_custom_nodes_centerpoints])
    
    prd = prd.merge(df_centerpoints[['LATITUDE', 'LONGITUDE']], left_on = 'FROM', right_index = True)
    prd = prd.merge(df_centerpoints[['LATITUDE', 'LONGITUDE']], left_on = 'TO', right_index = True, suffixes = ('_FROM', '_TO'))
    prd = prd.groupby(['TECHNOLOGY', 'YEAR', 'FROM', 'TO', 'LONGITUDE_FROM', 'LATITUDE_FROM', 'LONGITUDE_TO', 'LATITUDE_TO'],
                    as_index=False)['VALUE'].sum()
    
    if prd.empty:
        return
    
    scaler = MinMaxScaler()
    maxlinewidth = 3
    prd['line_width'] = (scaler.fit_transform(prd[['VALUE']]) * maxlinewidth).round(1)
    prd['line_width'].where(prd['line_width'] >= 0.3, 0.3, inplace = True)
    
    # get unique nodes to plot 
    nodes_to_plot = {}
    unique_nodes = set(prd["FROM"].to_list() + prd["TO"].to_list())
    for node in unique_nodes:
        nodes_to_plot[node] = [
            df_centerpoints.loc[node]["LATITUDE"],
            df_centerpoints.loc[node]["LONGITUDE"],
        ]
    # get plotting extent 
    lat_min = int(prd[['LATITUDE_FROM', 'LATITUDE_TO']].min().min() +-6)
    lat_max = int(prd[['LATITUDE_FROM', 'LATITUDE_TO']].max().max() +6)
    lon_min = int(prd[['LONGITUDE_FROM', 'LONGITUDE_TO']].min().min() +-6)
    lon_max = int(prd[['LONGITUDE_FROM', 'LONGITUDE_TO']].max().max() +6)
    extent = [lon_min, lon_max, lat_min, lat_max]
    
    # create new figure
    fig, ax = get_map(extent=extent)
    
    # generates all lines and map text
    df_year = prd.loc[prd['YEAR'] == int(year)]
    for y in df_year.index.unique():
        
        # get data to plot
        map_filter = df_year.loc[(df_year.index == y)]
        x1,y1 = map_filter['LONGITUDE_FROM'].iloc[0] , map_filter['LATITUDE_FROM'].iloc[0]
        x2,y2 = map_filter['LONGITUDE_TO'].iloc[0] , map_filter['LATITUDE_TO'].iloc[0]
        midco = midpoint(x1, y1, x2, y2)
        
        # add transmission line
        plot_map_trn_line(
            ax=ax,
            x1=x1,
            y1=y1,
            x2=x2,
            y2=y2,
            width=map_filter['line_width'].iloc[0]
        )
        
        # add capacity label
        plot_map_text(
            ax=ax,
            x=midco[0],
            y=midco[1],
            text=int(map_filter['VALUE'].iloc[0])
        )
        
    # add node labels 
    for _, node in enumerate(nodes_to_plot):
        plot_map_text(
            ax=ax,
            x=nodes_to_plot[node][1],
            y=nodes_to_plot[node][0],
            text=node
        )

    ax.set_title(f'Total Transmission Flow in {year} (PJ)', fontsize = '6')
    return fig.savefig(
        os.path.join(save_dir, f'TransmissionFlow{year}.jpg'), 
        dpi = 500, 
        bbox_inches="tight"
    )

if __name__ == '__main__':
        try:
            config_paths = ConfigPaths()
            config = ConfigFile('config')
            scenario_figs_dir = config_paths.scenario_figs_dir
            results_by_country = config.get('results_by_country')
            countries = config.get('geographic_scope')
            input_data = read_csv(config_paths.scenario_data_dir)
            result_data = read_csv(config_paths.scenario_results_dir)
            years = [config.get('endYear')]
            centerpoints = pd.read_csv(os.path.join(config_paths.input_data_dir, "centerpoints.csv"))
            custom_nodes = config.get('nodes_to_add')
            custom_nodes_centerpoints = pd.read_csv(os.path.join(config_paths.custom_nodes_dir, "centerpoints.csv"))
            
            main(input_data, 
                 result_data, 
                 centerpoints,
                 custom_nodes,
                 custom_nodes_centerpoints,
                 scenario_figs_dir, 
                 countries, 
                 results_by_country, 
                 years)
            
        except FileNotFoundError:
            print(f"Usage: python {sys.argv[0]} <input_data.csv> <result_data.csv> <scenario_figs_dir> <countries> <results_by_country> <years> <custom_nodes>")
            sys.exit(1)