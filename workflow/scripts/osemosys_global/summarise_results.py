import pandas as pd
import itertools
import os
import sys
import yaml
from OPG_configuration import ConfigFile, ConfigPaths
from visualisation import transform_ts, powerplant_filter
pd.set_option('mode.chained_assignment', None)


def main():
    '''Creates summaries of results'''

    config_paths = ConfigPaths()
    config = ConfigFile('config')

    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Check for output directory
    try:
        os.makedirs(scenario_result_summaries_dir)
    except FileExistsError:
        pass

    # SUMMARISE RESULTS
    headline_metrics()
    capacity_summary()
    generation_summary()
    trade_flows()


def renewables_filter(df):
    '''Function to filter and keep only renewable
    technologies'''
    renewables = ['BIO', 'CSP', 'GEO', 'HYD', 'SPV', 'WAS', 'WAV', 'WON', 'WOF']
    df = df[~df.TECHNOLOGY.str.contains('TRN')]
    df = df.loc[(df.TECHNOLOGY.str.startswith('PWR')) &
                (df.TECHNOLOGY.str[3:6].isin(renewables))
                ]
    return df


def fossil_filter(df):
    '''Function to filter and keep only fossil fuel
    technologies'''
    fossil_fuels = ['COA', 'COG', 'OCG', 'CCG', 'PET', 'OIL', 'OTH']
    df = df[~df.TECHNOLOGY.str.contains('TRN')]
    df = df.loc[(df.TECHNOLOGY.str.startswith('PWR')) &
                (df.TECHNOLOGY.str[3:6].isin(fossil_fuels))
                ]
    return df


def headline_metrics():
    '''Function to summarise metrics for:
       1. Total emissions
       2. Average RE share
       3. Average cost of electricity generation
       4. Total system cost
       5. Average fossil fuel share'''

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()

    # Fix path below to config_paths
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # GET RESULTS

    df_metrics = pd.DataFrame(columns=['Metric', 'Unit', 'Value'])
    df_metrics['Metric'] = ['Emissions',
                            'RE Share',
                            'Total System Cost',
                            'Cost of electricity',
                            'Fossil fuel share']
    df_metrics['Unit'] = ['Million tonnes of CO2-eq.',
                          '%',
                          'Billion $',
                          '$/MWh',
                          '%']

    # Emissions
    df_emissions = pd.read_csv(os.path.join(scenario_results_dir,
                                            'AnnualEmissions.csv'
                                            )
                               )
    emissions_total = df_emissions.VALUE.sum()
    df_metrics.loc[df_metrics['Metric'].str.startswith('Emissions'),
                   'Value'] = emissions_total.round(0)

    # Shares of RE and Fossil fuels
    df_shares = pd.read_csv(os.path.join(scenario_results_dir,
                                         'ProductionByTechnologyAnnual.csv'
                                         )
                            )
    df_shares = df_shares[~df_shares.TECHNOLOGY.str.contains('TRN')]
    gen_total = df_shares.loc[df_shares.TECHNOLOGY.str.startswith('PWR'),
                              'VALUE'].sum()

    # RE Share
    df_re_share = renewables_filter(df_shares)
    re_total = df_re_share.VALUE.sum()
    re_share = re_total / gen_total
    df_metrics.loc[df_metrics['Metric'].str.startswith('RE Share'),
                   'Value'] = (re_share*100).round(0)

    # Fossil Fuel Share
    df_ff_share = fossil_filter(df_shares)
    ff_total = df_ff_share.VALUE.sum()
    ff_share = ff_total / gen_total
    df_metrics.loc[df_metrics['Metric'].str.startswith('Fossil fuel share'),
                   'Value'] = (ff_share*100).round(0)

    # Total System Cost
    df_system_cost = pd.read_csv(os.path.join(scenario_results_dir,
                                              'TotalDiscountedCost.csv'
                                              )
                                 )
    system_cost_total = df_system_cost.VALUE.sum()
    df_metrics.loc[df_metrics['Metric'].str.startswith('Total System Cost'),
                   'Value'] = (system_cost_total/1000).round(0)

    # Cost of electricity generation
    df_demand = pd.read_csv(os.path.join(scenario_results_dir,
                                         'Demand.csv'
                                         )
                            )
    demand_total = df_demand.VALUE.sum()  # Total demand in TWh
    df_metrics.loc[df_metrics['Metric'].str.startswith('Cost of electricity'),
                   'Value'] = (system_cost_total/(demand_total*0.2778)).round(0)

    return df_metrics.to_csv(os.path.join(scenario_result_summaries_dir,
                                          'Metrics.csv'
                                          ),
                             index=None
                             )


def powerplant_summary(df):
    config_paths = ConfigPaths()
    input_data_dir = config_paths.input_data_dir
    name_colour_codes = pd.read_csv(os.path.join(input_data_dir,
                                                 'color_codes.csv'
                                                 ),
                                    encoding='latin-1')

    # Get colour mapping dictionary
    color_dict = dict([(i, n) for i, n
                       in zip(name_colour_codes.tech_id,
                              name_colour_codes.tech_name)])

    # Capacities
    df = df[~df.TECHNOLOGY.str.contains('TRN')]
    df = df[df.TECHNOLOGY.str.startswith('PWR')]
    df['NODE'] = (df['TECHNOLOGY'].str[6:9] +
                  '-' +
                  df['TECHNOLOGY'].str[9:11])
    df['LABEL'] = df['TECHNOLOGY'].str[3:6]
    df = df.replace(color_dict)
    
    return df


def capacity_summary():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Capacities
    df_capacities = pd.read_csv(os.path.join(scenario_results_dir,
                                             'TotalCapacityAnnual.csv'
                                             )
                                )

    df_capacities['NODE'] = (df_capacities['TECHNOLOGY'].str[6:9] +
                             '-' +
                             df_capacities['TECHNOLOGY'].str[9:11])
    df_capacities = powerplant_filter(df_capacities, country=None)
    df_capacities = df_capacities.groupby(['NODE', 'LABEL', 'YEAR'],
                                          as_index=False)['VALUE'].sum()
    df_capacities = df_capacities.sort_values(by=['YEAR',
                                                  'NODE',
                                                  'LABEL'])
    df_capacities['VALUE'] = df_capacities['VALUE'].round(2)

    return df_capacities.to_csv(os.path.join(scenario_result_summaries_dir,
                                             'Capacities.csv'
                                             ),
                                index=None
                                )


def generation_summary():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Capacities
    df_generation = pd.read_csv(os.path.join(scenario_results_dir,
                                             'ProductionByTechnology.csv'
                                             )
                                )
    df_generation = powerplant_filter(df_generation, country=None)
    df_generation = transform_ts(df_generation)
    df_generation = pd.melt(df_generation,
                            id_vars=['MONTH', 'HOUR', 'YEAR'],
                            value_vars=[x for x in df_generation.columns
                                        if x not in ['MONTH', 'HOUR', 'YEAR']],
                            value_name='VALUE')
    df_generation = df_generation[['YEAR', 'MONTH', 'HOUR', 'LABEL', 'VALUE']]

    return df_generation.to_csv(os.path.join(scenario_result_summaries_dir,
                                             'Generation.csv'
                                             ),
                                index=None
                                )


def trade_flows():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile('config')

    # Fix path below to config_paths
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir
    scenario_data_dir = config_paths.scenario_data_dir
    input_data_dir = config_paths.input_data_dir

    # GET TECHS TO PLOT

    df_gen = pd.read_csv(os.path.join(scenario_data_dir,
                                      'TECHNOLOGY.csv'))
    df_gen = df_gen[df_gen.VALUE.str.startswith('TRN')]
    interconnections = list(df_gen.VALUE.unique())

    # GET TIMESLICE DEFINITION

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

    model_start_year = config.get('startYear')
    model_end_year = config.get('endYear')
    years = list(range(model_start_year, model_end_year+1))

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

    df_ts_template = pd.DataFrame(list(itertools.product(interconnections,
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

    for each in dayparts_dict:
        df_ts_template.loc[(df_ts_template.HOUR > dayparts_dict[each][0]) &
                           (df_ts_template.HOUR <= dayparts_dict[each][1]),
                           'DAYPART'] = each

    # Trade flows
    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'TotalAnnualTechnologyActivityByMode.csv'
                                  )
                     )

    df['SEASON'] = df['TIMESLICE'].str[0:2]
    df['DAYPART'] = df['TIMESLICE'].str[2:]
    df['YEAR'] = df['YEAR'].astype(int)
    df.drop(['REGION', 'TIMESLICE'],
            axis=1,
            inplace=True)

    df = pd.merge(df,
                  df_ts_template,
                  how='left',
                  on=['TECHNOLOGY', 'SEASON', 'DAYPART', 'YEAR']).dropna()
    df['VALUE'] = (df['VALUE'].mul(1e6))/(df['DAYS'].mul(3600))

    df = df[['YEAR',
             'MONTH',
             'HOUR',
             'TECHNOLOGY',
             'MODE_OF_OPERATION',
             'VALUE']]
    df['MODE_OF_OPERATION'] = df['MODE_OF_OPERATION'].astype(int)
    df.loc[df['MODE_OF_OPERATION'] == 2, 'VALUE'] *= -1

    '''
    df['MODE_OF_OPERATION'].replace({1: 'NODE_1 to NODE_2',
                                     2: 'NODE_2 to NODE_1'},
                                    inplace=True)
    # Assign directions of trade flows

    df = df.pivot_table(index=['MONTH', 'HOUR', 'YEAR', 'TECHNOLOGY'],
                        columns='MODE_OF_OPERATION',
                        values='VALUE',
                        aggfunc='sum').reset_index().fillna(0)
    '''

    df['NODE_1'] = df.TECHNOLOGY.str[3:8]
    df['NODE_2'] = df.TECHNOLOGY.str[8:13]
    df.drop(columns=['TECHNOLOGY', 'MODE_OF_OPERATION'],
            axis=1,
            inplace=True)

    df['MONTH'] = pd.Categorical(df['MONTH'],
                                 categories=months,
                                 ordered=True)
    df = df.sort_values(by=['MONTH', 'HOUR'])
    df['VALUE'] = df['VALUE'].round(2)
    df = df[['YEAR',
             'MONTH',
             'HOUR',
             'NODE_1',
             'NODE_2',
             'VALUE']]

    return df.to_csv(os.path.join(scenario_result_summaries_dir,
                                  'TradeFlows.csv'
                                  ),
                     index=None
                     )



if __name__ == '__main__':
    main()
