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

#get paths from configuration file 
yaml_file = open("config.yaml")
parsed_yaml_file = yaml.load(yaml_file, Loader = yaml.FullLoader)
input_folder = os.path.join(parsed_yaml_file.get('outputDir'), 
                            parsed_yaml_file.get('scenario'), 
                            'results')
output_folder = os.path.join(parsed_yaml_file.get('outputDir'), 
                            parsed_yaml_file.get('scenario'), 
                            'figures')
model_folder = os.path.join(parsed_yaml_file.get('outputDir'), 
                            parsed_yaml_file.get('scenario'), 
                            'data')
data_folder = parsed_yaml_file.get('inputDir')


#input_folder = '../../osemosys_global_model/OsemosysGlobal/results'
#output_folder = '../../osemosys_global_model/OsemosysGlobal/figures'
#model_folder = '../../osemosys_global_model/OsemosysGlobal/data'
#data_folder = '../../data'

try:
    os.makedirs(output_folder)
except FileExistsError:
    pass

# Import tech-color codes table
name_color_codes = pd.read_csv(os.path.join(data_folder,
                                            'color_codes.csv'
                                            ),
                               encoding='latin-1')
tech_names = dict([(c, n) for c, n
                   in zip(name_color_codes.tech_id,
                          name_color_codes.tech_name)])
color_dict = dict([(n, c) for n, c
                   in zip(name_color_codes.tech_name,
                          name_color_codes.colour)])

# Create list of generation technologies
df_gen = pd.read_csv(os.path.join(model_folder,
                                  'TECHNOLOGY.csv'))
generation = list(df_gen.VALUE.unique())

# Import years
df_years = pd.read_csv(os.path.join(model_folder,
                                    'YEAR.csv'))
years = list(df_years.VALUE.unique())

# Import timeslice definition
seasons_months_days = pd.read_csv(os.path.join(data_folder,
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
dayparts_hours = pd.read_csv(os.path.join(data_folder,
                                          'ts_dayparts.csv'
                                          ),
                             encoding='latin-1')
dayparts_dict = dict(zip(dayparts_hours.daypart,
                         zip(dayparts_hours.start_hour,
                             dayparts_hours.end_hour))
                     )


def powerplant_filter(df):
    filtered_df = df[~df.TECHNOLOGY.str.contains('TRN')]
    filtered_df = filtered_df.loc[filtered_df.TECHNOLOGY.str[0:3] == 'PWR']
    filtered_df['TYPE'] = filtered_df.TECHNOLOGY.str[3:6]
    filtered_df['COUNTRY'] = filtered_df.TECHNOLOGY.str[6:9]
    filtered_df['LABEL'] = filtered_df['COUNTRY'] + '-' + filtered_df['TYPE']
    filtered_df.drop(['TECHNOLOGY', 'TYPE', 'COUNTRY'],
            axis=1,
            inplace=True)   
    return filtered_df


def transform_ts(df):
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
    # df = df.rename(columns = tech_names)
    return df


def plot_totalcapacity():
    df = pd.read_csv(os.path.join(input_folder,
                                  'TotalCapacityAnnual.csv'
                                  )
                     )
    df = powerplant_filter(df)
    df.VALUE = df.VALUE.astype('float64')
    df = df.groupby(['LABEL', 'YEAR'],
                    as_index=False)['VALUE'].sum()

    fig = px.bar(df,
                 x='YEAR',
                 y='VALUE',
                 color='LABEL',
                 # color_discrete_map=color_dict,
                 template='plotly_white',
                 labels={'YEAR': 'Year',
                         'VALUE': 'Gigawatts (GW)',
                         'LABEL': 'Country-Powerplant'})
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
    df = powerplant_filter(df)
    df.VALUE = df.VALUE.astype('float64')
    df = df.groupby(['LABEL', 'YEAR'],
                    as_index=False)['VALUE'].sum()

    fig = px.bar(df,
                 x='YEAR',
                 y='VALUE',
                 color='LABEL',
                 # color_discrete_map=color_dict,
                 template='plotly_white',
                 labels={'YEAR': 'Year',
                         'VALUE': 'Petajoules (PJ)',
                         'LABEL': 'Country-Powerplant'})
    fig.update_layout(
        font_family="Arial",
        font_size=14)
    fig.update_traces(marker_line_width=0,
                      opacity=0.8)

    return fig.write_html(os.path.join(output_folder,
                                       'GenerationAnnual.html'
                                       )
                          )


def plot_generation_hourly():
    df = pd.read_csv(os.path.join(input_folder,
                                  'ProductionByTechnology.csv'
                                  )
                     )
    df = powerplant_filter(df)
    df.VALUE = df.VALUE.astype('float64')
    df = transform_ts(df)

    fig = px.area(df,
                  x='HOUR',
                  y=[x for
                     x in df.columns
                     if x not in ['MONTH', 'HOUR']
                    ],
                  title='',
                  facet_col='MONTH',
                  facet_col_spacing=0.005,
                  # color_discrete_map=color_dict,
                  animation_frame='YEAR',
                  template='seaborn+plotly_white',
                  labels={
                      "variable": ""
                  }
                  )
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
    return fig.write_html(os.path.join(output_folder,
                                       'GenerationHourly.html'
                                       )
                          )


#plot_generation_hourly()
plot_totalcapacity()
#plot_generationannual()

