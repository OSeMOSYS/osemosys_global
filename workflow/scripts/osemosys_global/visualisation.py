import pandas as pd
pd.set_option('mode.chained_assignment', None)
import plotly as py
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import itertools
import os
import sys
import yaml
from OPG_configuration import ConfigFile, ConfigPaths

def main():
    '''Creates system level and country level graphs. '''

    config_paths = ConfigPaths()
    config = ConfigFile('config')
    scenario_figs_dir = config_paths.scenario_figs_dir
    results_by_country = config.get('results_by_country')

    # Check for output directory 
    try:
        os.makedirs(scenario_figs_dir)
    except FileExistsError:
        pass
    
    # Get system level results 
    plot_totalcapacity(country = None)
    plot_generationannual(country = None)
    plot_generation_hourly()

    # If producing by country results, check for folder structure 
    if results_by_country:
        countries = config.get('geographic_scope')
        for country in countries:
            try:
                os.makedirs(os.path.join(scenario_figs_dir, country))
            except FileExistsError:
                pass
    
            plot_totalcapacity(country = country)
            plot_generationannual(country = country)

def powerplant_filter(df, country = None):

    # CONFIGURATION PARAMETERS 
    config_paths = ConfigPaths()
    input_data_dir = config_paths.input_data_dir
    name_colour_codes = pd.read_csv(os.path.join(input_data_dir,
                                                'color_codes.csv'
                                                ),
                                   encoding='latin-1')

    # Get colour mapping dictionary
    color_dict = dict([(n, c) for n, c
                   in zip(name_colour_codes.tech_id,
                          name_colour_codes.colour)])

    filtered_df = df[~df.TECHNOLOGY.str.contains('TRN')]
    filtered_df = filtered_df.loc[filtered_df.TECHNOLOGY.str[0:3] == 'PWR']
    filtered_df['TYPE'] = filtered_df.TECHNOLOGY.str[3:6]
    filtered_df['COUNTRY'] = filtered_df.TECHNOLOGY.str[6:9]

    if country:
        filtered_df = filtered_df.loc[filtered_df['COUNTRY'] == country]
        filtered_df['LABEL'] = filtered_df['COUNTRY'] + '-' + filtered_df['TYPE']
    else:
        filtered_df['LABEL'] = filtered_df['TYPE']
    
    filtered_df['COLOR'] = filtered_df['TYPE'].map(color_dict)
    filtered_df.drop(['TECHNOLOGY', 'TYPE', 'COUNTRY'],
            axis=1,
            inplace=True)
    return filtered_df

def transform_ts(df):

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile('config')
    scenario_data_dir = config_paths.scenario_data_dir
    input_data_dir = config_paths.input_data_dir
    # years = range(config.get('startYear'), config.get('endYear') + 1)

    # GET TECHS TO PLOT

    df_gen = pd.read_csv(os.path.join(scenario_data_dir,
                                      'TECHNOLOGY.csv'))
    generation = list(df_gen.VALUE.unique())

    # GET TIMESLICE DEFENITION 
    '''
    seasons_months_days = pd.read_csv(os.path.join(input_data_dir,
                                               'ts_seasons.csv'
                                               ),
                                  encoding='latin-1')
    seasons_dict = dict([(m, s) for m, s
                         in zip(seasons_months_days.month_name,
                                seasons_months_days.season)])
    days_dict = dict([(m, d) for m, d
                      in zip(seasons_months_days.month_name,
                             seasons_months_days.days)])
    months = list(seasons_dict)
    hours = list(range(1, 25))

    dayparts_hours = pd.read_csv(os.path.join(input_data_dir,
                                              'ts_dayparts.csv'
                                              ),
                                 encoding='latin-1')
    dayparts_dict = dict(zip(dayparts_hours.daypart,
                             zip(dayparts_hours.start_hour,
                                 dayparts_hours.end_hour))
                         )
    '''

    seasons_raw = config.get('seasons')
    seasonsData = []

    for s, months in seasons_raw.items():
        for month in months:
            seasonsData.append([month, s]) 
    seasons_df = pd.DataFrame(seasonsData, 
                              columns=['month', 'season'])
    seasons_df = seasons_df.sort_values(by=['month']).reset_index(drop=True)
    dayparts_raw = config.get('dayparts')
    daypartData = []
    for dp, hr in dayparts_raw.items():
        daypartData.append([dp, hr[0], hr[1]])
    dayparts_df = pd.DataFrame(daypartData,
                               columns=['daypart', 'start_hour', 'end_hour'])
    timeshift = config.get('timeshift')
    dayparts_df['start_hour'] = dayparts_df['start_hour'].map(lambda x: apply_timeshift(x, timeshift))
    dayparts_df['end_hour'] = dayparts_df['end_hour'].map(lambda x: apply_timeshift(x, timeshift))

    month_names = {1: 'Jan',
                   2: 'Feb',
                   3: 'Mar',
                   4: 'Apr',
                   5: 'May',
                   6: 'Jun',
                   7: 'Jul',
                   8: 'Aug',
                   9: 'Sep',
                   10: 'Oct',
                   11: 'Nov',
                   12: 'Dec',
                   }

    days_per_month = {'Jan': 31,
                      'Feb': 28,
                      'Mar': 31,
                      'Apr': 30,
                      'May': 31,
                      'Jun': 30,
                      'Jul': 31,
                      'Aug': 31,
                      'Sep': 30,
                      'Oct': 31,
                      'Nov': 30,
                      'Dec': 31,
                      }

    seasons_df['month_name'] = seasons_df['month'].map(month_names)
    seasons_df['days'] = seasons_df['month_name'].map(days_per_month)
    seasons_df_grouped = seasons_df.groupby(['season'],
                                            as_index=False)['days'].sum()
    days_dict = dict(zip(list(seasons_df_grouped['season']),
                         list(seasons_df_grouped['days'])
                         )
                     )
    seasons_df['days'] = seasons_df['season'].map(days_dict)

    years = config.get_years()

    seasons_dict = dict(zip(list(seasons_df['month']),
                            list(seasons_df['season'])
                            )
                        )

    dayparts_dict = {i: [j, k]
                     for i, j, k
                     in zip(list(dayparts_df['daypart']),
                            list(dayparts_df['start_hour']),
                            list(dayparts_df['end_hour'])
                            )
                     }

    months = list(seasons_dict)
    hours = list(range(1, 25))

    # APPLY TRANSFORMATION

    df_ts_template = pd.DataFrame(list(itertools.product(generation,
                                                         months,
                                                         hours,
                                                         years)
                                       ),
                                  columns=['TECHNOLOGY',
                                           'MONTH',
                                           'HOUR',
                                           'YEAR']
                                  )

    df_ts_template = df_ts_template.sort_values(by=['TECHNOLOGY', 'YEAR'])
    df_ts_template['SEASON'] = df_ts_template['MONTH'].map(seasons_dict)
    df_ts_template['DAYS'] = df_ts_template['SEASON'].map(days_dict)
    df_ts_template['YEAR'] = df_ts_template['YEAR'].astype(int)
    df_ts_template = powerplant_filter(df_ts_template)

    for daypart in dayparts_dict:
        if dayparts_dict[daypart][0] > dayparts_dict[daypart][1]: # loops over 24hrs
            df_ts_template.loc[(df_ts_template['HOUR'] >= dayparts_dict[daypart][0]) |
                          (df_ts_template['HOUR'] < dayparts_dict[daypart][1]),
                          'DAYPART'] = daypart
        else:
            df_ts_template.loc[(df_ts_template['HOUR'] >= dayparts_dict[daypart][0]) &
                      (df_ts_template['HOUR'] < dayparts_dict[daypart][1]),
                      'DAYPART'] = daypart

    df['SEASON'] = df['TIMESLICE'].str[0:2]
    df['DAYPART'] = df['TIMESLICE'].str[2:]
    df['YEAR'] = df['YEAR'].astype(int)
    df.drop(['REGION', 'TIMESLICE'],
            axis=1,
            inplace=True)

    df = pd.merge(df,
                  df_ts_template,
                  how='left',
                  on=['LABEL', 'SEASON', 'DAYPART', 'YEAR']).dropna()
    df['VALUE'] = (df['VALUE'].mul(1e6))/(df['DAYS'].mul(3600))

    df = df.pivot_table(index=['MONTH', 'HOUR', 'YEAR'],
                        columns='LABEL',
                        values='VALUE',
                        aggfunc='sum').reset_index().fillna(0)

    df['MONTH'] = pd.Categorical(df['MONTH'],
                                 categories=months,
                                 ordered=True)
    df = df.sort_values(by=['MONTH', 'HOUR'])

    '''
    tech_names = dict([(c, n) for c, n
                   in zip(name_colour_codes.tech_id,
                          name_colour_codes.tech_name)])

    df = df.rename(columns = tech_names)
    '''
    return df


def plot_totalcapacity(country = None):

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_figs_dir = config_paths.scenario_figs_dir
    scenario_results_dir = config_paths.scenario_results_dir

    # GET RESULTS

    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'TotalCapacityAnnual.csv'
                                  )
                     )

    df = powerplant_filter(df, country)

    plot_colors = pd.Series(df['COLOR'].values, index=df['LABEL']).to_dict()
    df.VALUE = df.VALUE.astype('float64')
    df = df.groupby(['LABEL', 'YEAR'],
                    as_index=False)['VALUE'].sum()

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
    # fig.update_xaxes(type='category')
    fig.update_layout(
        font_family="Arial",
        font_size=14,
        legend_traceorder="reversed",
        title_x=0.5)
    fig['layout']['title']['font'] = dict(size=24)
    fig.update_traces(marker_line_width=0, opacity=0.8)

    if country:
        fig_file = os.path.join(scenario_figs_dir, country, 'TotalCapacityAnnual.html')
    else:
        fig_file = os.path.join(scenario_figs_dir, 'TotalCapacityAnnual.html')

    return fig.write_html(fig_file)

def plot_generationannual(country=None):

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_figs_dir = config_paths.scenario_figs_dir
    scenario_results_dir = config_paths.scenario_results_dir

    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'ProductionByTechnologyAnnual.csv'
                                  )
                     )

    df = powerplant_filter(df, country)

    plot_colors = pd.Series(df['COLOR'].values, index=df['LABEL']).to_dict()
    df.VALUE = df.VALUE.astype('float64')
    df = df.groupby(['LABEL', 'YEAR'],
                    as_index=False)['VALUE'].sum()

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
        fig_file = os.path.join(scenario_figs_dir, country, 'GenerationAnnual.html')
    else:
        fig_file = os.path.join(scenario_figs_dir, 'GenerationAnnual.html')

    return fig.write_html(fig_file)


def plot_generation_hourly():

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_figs_dir = config_paths.scenario_figs_dir
    scenario_results_dir = config_paths.scenario_results_dir

    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'ProductionByTechnology.csv'
                                  )
                     )
    df = powerplant_filter(df, country=False)
    plot_colors = pd.Series(df['COLOR'].values, index=df['LABEL']).to_dict()
    df.VALUE = df.VALUE.astype('float64')
    df = transform_ts(df)

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
        yaxis_title = 'Petajoules (PJ)')
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
    return fig.write_html(os.path.join(scenario_figs_dir,
                                       'GenerationHourly.html'
                                       )
                          )

def apply_timeshift(x, timeshift):
    '''Applies timeshift to organize dayparts.
    
    Args:
        x = Value between 0-24
        timeshift = value offset from UTC (-11 -> +12)'''

    x += timeshift
    if x > 23:
        return x - 24
    elif x < 0:
        return x + 24
    else:
        return x

if __name__ == '__main__':
    main()
