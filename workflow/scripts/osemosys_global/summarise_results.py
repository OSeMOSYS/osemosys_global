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

    # # SUMMARISE RESULTS
    headline_metrics()
    capacity_summary()
    generation_summary()
    generation_by_node_summary()
    trade_flows()
    
    # UPDATED METRICS
    system_cost_by_node()
    new_capacity_summary()
    new_capacity_summary_trn()
    investment_summary()
    investment_summary_trn()
    marginal_costs()


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
                            'Renewable energy share',
                            'Total system cost',
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
    df_metrics.loc[df_metrics['Metric'].str.startswith('Renewable energy share'),
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
    df_metrics.loc[df_metrics['Metric'].str.startswith('Total system cost'),
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
    
def new_capacity_summary():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Capacities
    df_capacities = pd.read_csv(os.path.join(scenario_results_dir,
                                             'NewCapacity.csv'
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
                                             'New_Capacities_Powerplants.csv'
                                             ),
                                index=None
                                )
    
def new_capacity_summary_trn():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Capacities
    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'NewCapacity.csv'
                                  )
                     )
    df = df[df.TECHNOLOGY.str.startswith('TRN')]
    interconnections = list(df.TECHNOLOGY.unique())
    if len(interconnections) > 0:
        df['NODE_1'] = df.TECHNOLOGY.str[3:8]
        df['NODE_2'] = df.TECHNOLOGY.str[8:13]
        
        min_year = df['YEAR'].min()
        df = df[df['YEAR'] > min_year]
    
        df = df[['YEAR',
                 'NODE_1',
                 'NODE_2',
                 'VALUE']]
    else:
        df = pd.DataFrame(columns=['YEAR',
                                   'NODE_1',
                                   'NODE_2',
                                   'VALUE']
                         )

    return df.to_csv(os.path.join(scenario_result_summaries_dir,
                                  'New_Capacities_Interconnectors.csv'
                                  ),
                     index=None
                     )


def investment_summary():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Capacities
    df_capacities = pd.read_csv(os.path.join(scenario_results_dir,
                                             'CapitalInvestment.csv'
                                             )
                                )

    df_capacities['NODE'] = (df_capacities['TECHNOLOGY'].str[6:9] +
                             '-' +
                             df_capacities['TECHNOLOGY'].str[9:11])
    df_capacities = powerplant_filter(df_capacities, country=None)
    df_capacities = df_capacities.groupby(['NODE', 'LABEL', 'YEAR'],
                                          as_index=False)['VALUE'].sum()
    min_year = df_capacities['YEAR'].min()
    df_capacities = df_capacities[df_capacities['YEAR'] > min_year]
    df_capacities = df_capacities.sort_values(by=['YEAR',
                                                  'NODE',
                                                  'LABEL'])
    df_capacities['VALUE'] = df_capacities['VALUE'].round(2)

    return df_capacities.to_csv(os.path.join(scenario_result_summaries_dir,
                                             'Investment_Summary_Powerplants.csv'
                                             ),
                                index=None
                                )
    
def investment_summary_trn():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Capacities
    df = pd.read_csv(os.path.join(scenario_results_dir,
                                  'CapitalInvestment.csv'
                                  )
                     )
    df = df[df.TECHNOLOGY.str.startswith('TRN')]
    interconnections = list(df.TECHNOLOGY.unique())
    if len(interconnections) > 0:
        df['NODE_1'] = df.TECHNOLOGY.str[3:8]
        df['NODE_2'] = df.TECHNOLOGY.str[8:13]
        
        min_year = df['YEAR'].min()
        df = df[df['YEAR'] > min_year]
    
        df = df[['YEAR',
                 'NODE_1',
                 'NODE_2',
                 'VALUE']]
    else:
        df = pd.DataFrame(columns=['YEAR',
                                   'NODE_1',
                                   'NODE_2',
                                   'VALUE']
                         )

    return df.to_csv(os.path.join(scenario_result_summaries_dir,
                                  'Investment_Summary_Interconnectors.csv'
                                  ),
                     index=None
                     )


def generation_summary():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir

    # Generation
    df_generation = pd.read_csv(os.path.join(scenario_results_dir,
                                             'ProductionByTechnology.csv'
                                             )
                                )
    df_generation = powerplant_filter(df_generation, country=None)
    df_generation = df_generation.loc[df_generation['FUEL'].str.startswith('ELC')]
    df_generation = transform_ts(df_generation)
    df_generation = pd.melt(df_generation,
                            id_vars=['MONTH', 'HOUR', 'YEAR'],
                            value_vars=[x for x in df_generation.columns
                                        if x not in ['MONTH', 'HOUR', 'YEAR']],
                            value_name='VALUE')
    df_generation['VALUE'] = df_generation['VALUE'].round(2)
    df_generation = df_generation[['YEAR', 'MONTH', 'HOUR', 'LABEL', 'VALUE']]

    return df_generation.to_csv(os.path.join(scenario_result_summaries_dir,
                                             'Generation.csv'
                                             ),
                                index=None
                                )


def generation_by_node_summary():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile('config')
    scenario_results_dir = config_paths.scenario_results_dir
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir
    scenario_data_dir = config_paths.scenario_data_dir
    input_data_dir = config_paths.input_data_dir

    # Generation
    df_gen_by_node = pd.read_csv(os.path.join(scenario_results_dir,
                                             'ProductionByTechnology.csv'
                                             )
                                )
    # GET TECHS TO PLOT
    generation = list(df_gen_by_node.TECHNOLOGY.unique())
    df_gen_by_node['NODE'] = (df_gen_by_node['TECHNOLOGY'].str[6:11])
    df_gen_by_node = powerplant_filter(df_gen_by_node, country=None)
    df_gen_by_node = df_gen_by_node.loc[df_gen_by_node['FUEL'].str.startswith('ELC')]
    # df_generation = transform_ts(df_generation)

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

    hours_dict = {i: abs(k-j)
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

    df_ts_template = df_ts_template.drop_duplicates()

    df_gen_by_node['SEASON'] = df_gen_by_node['TIMESLICE'].str[0:2]
    df_gen_by_node['DAYPART'] = df_gen_by_node['TIMESLICE'].str[2:]
    df_gen_by_node['YEAR'] = df_gen_by_node['YEAR'].astype(int)
    df_gen_by_node.drop(['REGION', 'TIMESLICE'],
                        axis=1,
                        inplace=True)

    df_gen_by_node = pd.merge(df_gen_by_node,
                              df_ts_template,
                              how='left',
                              on=['LABEL', 'SEASON', 'DAYPART', 'YEAR']).dropna()
    df_gen_by_node['HOUR_COUNT'] = df_gen_by_node['DAYPART'].map(hours_dict)
    df_gen_by_node['VALUE'] = (df_gen_by_node['VALUE'].mul(1e6))/(df_gen_by_node['DAYS']*df_gen_by_node['HOUR_COUNT'].mul(3600))

    df_gen_by_node = df_gen_by_node.pivot_table(index=['MONTH', 'HOUR', 'YEAR', 'NODE'],
                                              columns='LABEL',
                                              values='VALUE',
                                              aggfunc='sum').reset_index().fillna(0)

    df_gen_by_node['MONTH'] = pd.Categorical(df_gen_by_node['MONTH'],
                                            categories=months,
                                            ordered=True)
    df_gen_by_node = df_gen_by_node.sort_values(by=['MONTH', 'HOUR'])
    '''
    df_generation = pd.melt(df_generation,
                            id_vars=['MONTH', 'HOUR', 'YEAR', 'NODE'],
                            value_vars=[x for x in df_generation.columns
                                        if x not in ['MONTH', 'HOUR', 'YEAR', 'NODE']],
                            value_name='VALUE')
    '''
    cols_round = [x for x in df_gen_by_node.columns
                  if x not in ['MONTH', 'HOUR', 'YEAR', 'NODE']]
    df_gen_by_node[cols_round] = df_gen_by_node[cols_round].round(2)
    # df_generation = df_generation[['YEAR', 'MONTH', 'HOUR', 'NODE', 'LABEL', 'VALUE']]

    return df_gen_by_node.to_csv(os.path.join(scenario_result_summaries_dir,
                                             'Generation_By_Node.csv'
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

    if len(interconnections) > 0:
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

        hours_dict = {i: abs(k-j)
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

        for daypart in dayparts_dict:
            if dayparts_dict[daypart][0] > dayparts_dict[daypart][1]: # loops over 24hrs
                df_ts_template.loc[(df_ts_template['HOUR'] >= dayparts_dict[daypart][0]) |
                            (df_ts_template['HOUR'] < dayparts_dict[daypart][1]),
                            'DAYPART'] = daypart
            else:
                df_ts_template.loc[(df_ts_template['HOUR'] >= dayparts_dict[daypart][0]) &
                        (df_ts_template['HOUR'] < dayparts_dict[daypart][1]),
                        'DAYPART'] = daypart

        df_ts_template = df_ts_template.drop_duplicates()

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

        df['HOUR_COUNT'] = df['DAYPART'].map(hours_dict)
        df['VALUE'] = (df['VALUE'].mul(1e6))/(df['DAYS']*df['HOUR_COUNT'].mul(3600))

        df = df[['YEAR',
                'MONTH',
                'HOUR',
                'TECHNOLOGY',
                'MODE_OF_OPERATION',
                'VALUE']]
        df['MODE_OF_OPERATION'] = df['MODE_OF_OPERATION'].astype(int)
        df.loc[df['MODE_OF_OPERATION'] == 2, 'VALUE'] *= -1

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
    else:
        df = pd.DataFrame(columns=['YEAR',
                                   'MONTH',
                                   'HOUR',
                                   'NODE_1',
                                   'NODE_2',
                                   'VALUE']
                         )

    return df.to_csv(os.path.join(scenario_result_summaries_dir,
                                  'TradeFlows.csv'
                                  ),
                     index=None
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

'''
    + Total system costs for each node [$]
    + Existing, new and decommissioned capacity in electricity generation, 
    storage, and transmission technology for each year and node [GW (and GWh in 
    the case of storage)]
    + Electricity generation, charge and discharge and trade by technology for 
    each year and node [GWh]
    + Emissions of each fossil generation technology for each year and node [tCO2]
    x Implicit electricity price for each year and node [$/MWh]
    x Implicit carbon price [$/t]
    x Share in electricity production and generation capacity [%]
    + Investment volumes by technology [$]
    x Energetic storage losses node and technology [GWh]
    - Capacity factors by node and technology [%]
    - Curtailment by node and technology [GWh]
    + Fossil fuel consumption [t for coal, bcf for gas]
    x Early decommissioning of fossil generation [GW]
'''

def system_cost_by_node():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile('config')
    scenario_results_dir = config_paths.scenario_results_dir
    #scenario_results_dir = '/Users/adminuser/Documents/repositories/feo-esmod-osemosys/workflow/scripts/osemosys_global/../../../results/Indonesia_BA/results'
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir
    #scenario_result_summaries_dir = '/Users/adminuser/Documents/repositories/feo-esmod-osemosys/workflow/scripts/osemosys_global/../../../results/Indonesia_BA/result_summaries'
    scenario_data_dir = config_paths.scenario_data_dir
    #scenario_data_dir = '/Users/adminuser/Documents/repositories/feo-esmod-osemosys/workflow/scripts/osemosys_global/../../../results/Indonesia_BA/data'
    input_data_dir = config_paths.input_data_dir
    
    penalty = config.get('emission_penalty')
    
    '''
    df = pd.read_csv(os.path.join(scenario_data_dir,
                                  'TECHNOLOGY.csv'))
    df.rename(columns = {'VALUE': 'TECHNOLOGY'},
              inplace=True)
    
    df['VALUE'] = 0
    df.set_index('TECHNOLOGY', inplace=True)
    '''
    df = pd.DataFrame(columns=['TECHNOLOGY','YEAR','VALUE'])
    
    # System costs by node
    
    # Investment costs
    df_inv = pd.read_csv(os.path.join(scenario_results_dir,
                                      'CapitalInvestment.csv'
                                      ))
    df_inv = df_inv[~(df_inv['TECHNOLOGY'].str.startswith('MIN')) &
                    ~(df_inv['TECHNOLOGY'].str.startswith('RNW')) &
                    ~(df_inv['TECHNOLOGY'].str.startswith('TRN'))]
    
    # Storage costs
    if os.path.exists(os.path.join(scenario_results_dir,
                                   'NewStorageCapacity.csv'
                                   )):
        df_sto = pd.read_csv(os.path.join(scenario_results_dir,
                                          'NewStorageCapacity.csv'
                                          ))
        df_sto_cost = pd.read_csv(os.path.join(scenario_data_dir,
                                               'CapitalCostStorage.csv'
                                               ))
        df_sto_cost.rename(columns={'VALUE':'COST'},
                        inplace=True)
        
        df_sto = pd.merge(df_sto, df_sto_cost,
                        on=['REGION','STORAGE', 'YEAR'],
                        how='left')
        df_sto['VALUE'] = df_sto['VALUE'] * df_sto_cost['COST']
        df_sto['STORAGE'] = 'PWR' + df_sto['STORAGE']
        df_sto.rename(columns={'STORAGE':'TECHNOLOGY'},
                    inplace=True)
        df_sto = df_sto[['REGION','TECHNOLOGY','YEAR','VALUE']]
    else:
        df_sto = pd.DataFrame(columns=['REGION','TECHNOLOGY','YEAR','VALUE'])
    
    # Fixed O&M
    df_fom = pd.read_csv(os.path.join(scenario_results_dir,
                                      'AnnualFixedOperatingCost.csv'
                                      ))
    
    # Variable O&M
    df_vom = pd.read_csv(os.path.join(scenario_results_dir,
                                      'AnnualVariableOperatingCost.csv'
                                      ))
    df_vom = df_vom[~(df_vom['TECHNOLOGY'].str.startswith('MIN'))]
    
    # Emissions penalty
    df_emi = pd.read_csv(os.path.join(scenario_results_dir,
                                      'AnnualTechnologyEmission.csv'
                                      ))
    df_pen = pd.read_csv(os.path.join(scenario_data_dir,
                                      'EmissionsPenalty.csv'
                                      ))
    df_pen.rename(columns={'VALUE': 'PENALTY'},
                  inplace=True)
    
    df_emi_pen = pd.merge(df_emi, df_pen,
                          how='left',
                          on=['REGION', 'EMISSION', 'YEAR'])
    df_emi_pen['VALUE'] = df_emi_pen['VALUE'] * df_emi_pen['PENALTY']
    df_emi_pen = df_emi_pen.groupby(['REGION', 'TECHNOLOGY', 'YEAR'],
                                    as_index=False)['VALUE'].sum()
    df_emi_pen.dropna(inplace=True)
    
    for each_df in [df_inv, df_fom, df_vom, df_sto, df_emi_pen]:
        each_df = each_df[['TECHNOLOGY',
                           'YEAR',
                           'VALUE']].fillna(0)
        #each_df.set_index(['TECHNOLOGY', 'YEAR'], inplace=True)
        #df = df.add(each_df, fill_value=0)
        
        df = pd.concat([df, each_df])

    df = df.groupby(['TECHNOLOGY','YEAR'],
                    as_index=False)['VALUE'].sum()
    
    df.reset_index(inplace=True)
    df = df[~(df['TECHNOLOGY'].str.startswith('MIN')) & 
            ~(df['TECHNOLOGY'].str.startswith('RNW')) &
            ~(df['TECHNOLOGY'].str.startswith('TRN'))]
    df['NODE'] = df['TECHNOLOGY'].str[6:11]
    df = df.groupby(['NODE',
                     'YEAR'],
                    as_index=False)['VALUE'].sum()
    #df = df[['NODE',
    #         'VALUE']]
    
    # Summarise UseByTechnologyAnnual for all powerplants
    df_use = pd.read_csv(os.path.join(scenario_results_dir,
                                      'TotalAnnualTechnologyActivityByMode.csv'
                                      ))
    df_use = df_use.groupby(['TECHNOLOGY', 'MODE_OF_OPERATION', 'YEAR'],
                            as_index=False)['VALUE'].sum()
    df_use = df_use.loc[(df_use['TECHNOLOGY'].str.startswith('PWR')) & 
                        ~(df_use['TECHNOLOGY'].str.startswith('PWRBAT')) &
                        ~(df_use['TECHNOLOGY'].str.startswith('PWRTRN'))]
    df_use['VALUE'] = df_use['VALUE'].round(4)
    
    # Get InputActivityRatios to calculate use by mode of operation
    df_iar = pd.read_csv(os.path.join(scenario_data_dir,
                                      'InputActivityRatio.csv'
                                      ))
    df_iar = df_iar[['TECHNOLOGY','FUEL','MODE_OF_OPERATION','YEAR','VALUE']]
    df_iar.rename(columns={'VALUE':'IAR'},
                  inplace=True)
    
    # Calculate use from activity and IAR for each technology
    df_use = pd.merge(df_use, df_iar,
                      on=['TECHNOLOGY', 'MODE_OF_OPERATION', 'YEAR'],
                      how='left')
    df_use['USE'] = df_use['VALUE']*df_use['IAR'] 
    
    # Get OutputActivityRatio to get fuel costs
    df_oar = pd.read_csv(os.path.join(scenario_data_dir,
                                      'OutputActivityRatio.csv'
                                      ))
    df_oar = df_oar[['TECHNOLOGY',
                     'FUEL',
                     'MODE_OF_OPERATION',
                     'YEAR']]
    
    # Get VAR for each technology
    df_var = pd.read_csv(os.path.join(scenario_data_dir,
                                      'VariableCost.csv'
                                      ))
    df_var = df_var[['TECHNOLOGY',
                     'MODE_OF_OPERATION',
                     'YEAR',
                     'VALUE']]
    df_var['VALUE'] = df_var['VALUE'].round(2)
    df_var.rename(columns={'VALUE':'VAR'},
                  inplace=True)
    
    # Get EAR for each technology
    df_ear = pd.read_csv(os.path.join(scenario_data_dir,
                                      'EmissionActivityRatio.csv'
                                      ))
    df_ear = df_ear[['TECHNOLOGY',
                     'MODE_OF_OPERATION',
                     'YEAR',
                     'VALUE']]
    df_ear['VALUE'] = df_ear['VALUE'].round(4)
    df_ear.rename(columns={'VALUE':'EAR'},
                  inplace=True)
    
    # Combine VAR and EAR for each technology
    df_var_ear = pd.merge(df_var, df_ear,
                          on=['TECHNOLOGY', 'MODE_OF_OPERATION', 'YEAR'],
                          how='outer')
    
    # Get fuel for each MIN technology from OAR 
    df_var_ear_oar = pd.merge(df_oar, df_var_ear,
                              on=['TECHNOLOGY', 'MODE_OF_OPERATION', 'YEAR'],
                              how='left')
    df_var_ear_oar = df_var_ear_oar[['FUEL', 'YEAR', 'VAR', 'EAR']].drop_duplicates().dropna()
    
    # Consolidate all data in single table (VAR, OAR, EAR, USE)
    df_all = pd.merge(df_use, df_var_ear_oar,
                      on=['FUEL', 'YEAR'],
                      how='outer')
    
    df_all['FUEL_COST'] = df_all['USE']*df_all['VAR'] 
    #df_all['EMISSIONS'] = df_all['USE']*df_all['EAR']
    df_all.dropna(inplace=True) 
    df_all = df_all[df_all['TECHNOLOGY'].str.startswith('PWR')]
    df_all['LABEL'] = df_all['TECHNOLOGY'].str[3:6]
    df_all['NODE'] = df_all['TECHNOLOGY'].str[6:11]
        
    # Calculate summary tables
    df_fuel_cost = df_all.groupby(['NODE','LABEL','YEAR'],
                                  as_index=False)['FUEL_COST'].sum()
    df_fuel_cost = df_fuel_cost[df_fuel_cost['FUEL_COST'] > 0]
    
    '''
    df_emissions = df_all.groupby(['NODE','LABEL','YEAR'],
                                  as_index=False)['EMISSIONS'].sum()
    df_emissions = df_emissions[~(df_emissions['LABEL'].str.startswith('CCS'))]
    df_emissions = df_emissions[df_emissions['EMISSIONS'] > 0]
    df_emissions.to_csv(os.path.join(scenario_result_summaries_dir,
                                     'AnnualEmissionsByNode.csv'
                                     ),
                        index=None
                        )
    df_emissions['PENALTY'] = df_emissions['EMISSIONS'] * penalty
    '''
    
    df_fuel_use = df_all.groupby(['NODE','LABEL'],
                                  as_index=False)['USE'].sum()
    df_fuel_use = df_fuel_use[df_fuel_use['LABEL'].isin(['COA',
                                                         'CCG',
                                                         'OCG'])]
    df_fuel_use['LABEL'] = df_fuel_use['LABEL'].replace(['CCG','OCG'], 'GAS')
    df_fuel_use = df_fuel_use.groupby(['NODE','LABEL'],
                                      as_index=False)['USE'].sum()
    df_fuel_use = df_fuel_use[df_fuel_use['USE'] > 0]
    df_fuel_use = df_fuel_use.pivot(index='NODE',
                                    columns='LABEL',
                                    values='USE').reset_index().fillna(0)
    if 'COA' in df_fuel_use.columns:
        df_fuel_use['COA'] = df_fuel_use['COA'] / 19 # Energy content of 19 MJ/kg
    if 'GAS' in df_fuel_use.columns:
        df_fuel_use['GAS'] = df_fuel_use['GAS'] * 0.9478 # PJ to bcf of Natural Gas
    df_fuel_use.to_csv(os.path.join(scenario_result_summaries_dir,
                                    'FuelUse.csv'
                                    ),
                       index=None
                       )
    
    df_fuel_cost = df_fuel_cost.groupby(['NODE','YEAR'],
                                        as_index=False)['FUEL_COST'].sum()
    df = pd.merge(df, df_fuel_cost,
                  on=['NODE', 'YEAR'],
                  how='outer')
    
    #df_emi_pen = df_emissions.groupby(['NODE','YEAR'],
    #                                  as_index=False)['PENALTY'].sum()
    #df = pd.merge(df, df_emi_pen,
    #              on=['NODE', 'YEAR'],
    #              how='outer')
    df.fillna(0,
              inplace=True)
    #df['SYSTEM_COST'] = df['VALUE'] + df['FUEL_COST'] + df['PENALTY']
    df['SYSTEM_COST'] = df['VALUE'] + df['FUEL_COST']
    df = df.groupby(['NODE'],
                    as_index=False)['SYSTEM_COST'].sum()
        
        
    return df.to_csv(os.path.join(scenario_result_summaries_dir,
                                  'SystemCostByNode.csv'
                                  ),
                     index=None
                     )


def marginal_costs():
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile('config')
    scenario_results_dir = config_paths.scenario_results_dir
    scenario = config.get('scenario')
    scenario_dir = config_paths.scenario_dir
    #scenario = 'ASEAN_v4_APG_LC'
    #scenario_results_dir = '/Users/adminuser/Documents/repositories/feo-esmod-osemosys/results/' + scenario
    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir
    
    
    duals = []
    
    with open(os.path.join(scenario_dir,
                           scenario + '.attr')) as sol_file:
        for line in sol_file:
            if line.startswith('EBa11'):
                ts = line.split(' ')[0].split(',')[1]
                fuel = line.split(' ')[0].split(',')[2]
                year = int(line.split(' ')[0].split(',')[3].split(')')[0])
                value = float(line.split(' ')[1])
                if fuel.startswith('ELC'):
                    if fuel.endswith('02'):
                        duals.append([ts, fuel, year, round(value, 2)])
    
    df_duals = pd.DataFrame(duals,
                            columns=['TS',
                                     'FUEL',
                                     'YEAR',
                                     'VALUE'])
    df_duals['SEASON'] = df_duals['TS'].str[:2]
    df_duals['DAYPART'] = df_duals['TS'].str[2:]
    months = list(range(1, 13))
    hours = list(range(1, 25))
    
    # Create DataFrame scaffold for dual values
    df_duals_final = pd.DataFrame(list(itertools.product(df_duals['FUEL'].unique(),
                                                         months,
                                                         hours,
                                                         df_duals['YEAR'].unique())
                                       ),
                                       columns=['FUEL',
                                                'MONTH',
                                                'HOUR',
                                                'YEAR']
                                  )
    # Create dictionaries of seasons and dayparts
    seasons_raw = config.get('seasons')
    seasons_dict = {}

    for s, months in seasons_raw.items():
        for month in months:
            seasons_dict[month] = s
    
    dayparts_raw = config.get('dayparts')
    dayparts_dict = {}
    for dp, hours in dayparts_raw.items():
        for hour in range(hours[0], hours[1]):
            dayparts_dict[hour+1] = dp

    # Create SEASON and DAYPART columns for each hour and month
    df_duals_final['SEASON'] = df_duals_final['MONTH'].map(seasons_dict)
    df_duals_final['DAYPART'] = df_duals_final['HOUR'].map(dayparts_dict)
    
    df_duals_final = pd.merge(df_duals_final, df_duals,
                              how='left',
                              on=['FUEL', 'SEASON', 'DAYPART', 'YEAR'])
    
    df_duals_final['NODE'] = df_duals_final['FUEL'].str[3:8]
    
    # Filter columns for final DataFrame
    df_duals_final = df_duals_final[['NODE',
                                     'MONTH',
                                     'HOUR',
                                     'YEAR',
                                     'VALUE']]
    
    # Convert $mn/PJ to $/MWh i.e. 3.6
    df_duals_final['VALUE'] = df_duals_final['VALUE'].mul(3.6)
    
    # Apply timeshift
    timeshift = config.get('timeshift')
    #df_duals_final['HOUR'] = df_duals_final['HOUR'].map(lambda x: apply_timeshift(x, timeshift))    
    
    #print(df_duals_final)

    return df_duals_final.to_csv(os.path.join(scenario_result_summaries_dir,
                                              'SRMC.csv'),
                                 index=None)


if __name__ == '__main__':
    main()
