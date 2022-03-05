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

# PATH CONSTANTS 

_YAML_FILE = open(os.path.join(os.path.dirname(__file__), '../../..',
                              'config/config.yaml'))
_PARSED_YAML_FILE = yaml.load(_YAML_FILE, Loader = yaml.FullLoader)

_INPUT_FOLDER = os.path.join(os.path.dirname(__file__), '../../..',
                            _PARSED_YAML_FILE.get('outputDir'), 
                            _PARSED_YAML_FILE.get('scenario'), 
                            'results')
_OUTPUT_FOLDER = os.path.join(os.path.dirname(__file__), '../../..',
                            _PARSED_YAML_FILE.get('outputDir'), 
                            _PARSED_YAML_FILE.get('scenario'), 
                            'figures')
_MODEL_FOLDER = os.path.join(os.path.dirname(__file__), '../../..',
                            _PARSED_YAML_FILE.get('outputDir'), 
                            _PARSED_YAML_FILE.get('scenario'), 
                            'data')
_DATA_FOLDER = os.path.join(os.path.dirname(__file__), '../../..',
                           _PARSED_YAML_FILE.get('inputDir'), 'data')

# IMPORT TECH COLOUR CODES

_NAME_COLOUR_CODES = pd.read_csv(os.path.join(_DATA_FOLDER,
                                                'color_codes.csv'
                                                ),
                                   encoding='latin-1')

def main():
    '''Creates system level and country level graphs. '''

    # Check for output directory 
    try:
        os.makedirs(_OUTPUT_FOLDER)
    except FileExistsError:
        pass
    
    # Get system level results 
    plot_generation_hourly()
    plot_totalcapacity(country = None)
    plot_generationannual(country = None)

    # Flag if to produce results per country 
    results_by_country = _PARSED_YAML_FILE.get('results_by_country')

    # If producing by country results, check for folder structure 
    if results_by_country:
        countries = _PARSED_YAML_FILE.get('geographic_scope')
        for country in countries:
            try:
                os.makedirs(os.path.join(_OUTPUT_FOLDER, country))
            except FileExistsError:
                pass
    
            plot_totalcapacity(country = country)
            plot_generationannual(country = country)

def powerplant_filter(df, country = None):

    # Get colour mapping dictionary
    color_dict = dict([(n, c) for n, c
                   in zip(_NAME_COLOUR_CODES.tech_id,
                          _NAME_COLOUR_CODES.colour)])

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

    # GET GENERATION TECHS

    df_gen = pd.read_csv(os.path.join(_MODEL_FOLDER,
                                      'TECHNOLOGY.csv'))
    generation = list(df_gen.VALUE.unique())

    # GET YEARS 

    years = range(
        _PARSED_YAML_FILE.get('startYear'),
        _PARSED_YAML_FILE.get('endYear') + 1,
    )

    # GET TIMESLICE DEFENITION 

    seasons_months_days = pd.read_csv(os.path.join(_DATA_FOLDER,
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

    dayparts_hours = pd.read_csv(os.path.join(_DATA_FOLDER,
                                              'ts_dayparts.csv'
                                              ),
                                 encoding='latin-1')
    dayparts_dict = dict(zip(dayparts_hours.daypart,
                             zip(dayparts_hours.start_hour,
                                 dayparts_hours.end_hour))
                         )

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
    df_ts_template['DAYS'] = df_ts_template['MONTH'].map(days_dict)
    df_ts_template['SEASON'] = df_ts_template['MONTH'].map(seasons_dict)
    df_ts_template['YEAR'] = df_ts_template['YEAR'].astype(int)
    df_ts_template = powerplant_filter(df_ts_template)

    for each in dayparts_dict:
        df_ts_template.loc[(df_ts_template.HOUR > dayparts_dict[each][0]) &
                           (df_ts_template.HOUR <= dayparts_dict[each][1]),
                           'DAYPART'] = each

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
                   in zip(_NAME_COLOUR_CODES.tech_id,
                          _NAME_COLOUR_CODES.tech_name)])

    df = df.rename(columns = tech_names)
    '''
    return df


def plot_totalcapacity(country = None):
    df = pd.read_csv(os.path.join(_INPUT_FOLDER,
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
        fig_file = os.path.join(_OUTPUT_FOLDER, country, 'TotalCapacityAnnual.html')
    else:
        fig_file = os.path.join(_OUTPUT_FOLDER, 'TotalCapacityAnnual.html')

    return fig.write_html(fig_file)

def plot_generationannual(country=None):
    df = pd.read_csv(os.path.join(_INPUT_FOLDER,
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
        fig_file = os.path.join(_OUTPUT_FOLDER, country, 'GenerationAnnual.html')
    else:
        fig_file = os.path.join(_OUTPUT_FOLDER, 'GenerationAnnual.html')

    return fig.write_html(fig_file)


def plot_generation_hourly():
    df = pd.read_csv(os.path.join(_INPUT_FOLDER,
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
                      "variable": ""})
    fig.update_layout(
        legend_traceorder="reversed",
        title_x=0.5)
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
    return fig.write_html(os.path.join(_OUTPUT_FOLDER,
                                       'GenerationHourly.html'
                                       )
                          )

if __name__ == '__main__':
    main()
