import pandas as pd
import itertools
import os
import sys
import yaml
from OPG_configuration import ConfigFile, ConfigPaths
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
    # scenario_results_dir = '/home/abhi/osemosys_global/results/osemosys-global/results'
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
    df['POWERPLANT'] = df['TECHNOLOGY'].str[3:6]
    df = df.groupby(['NODE', 'POWERPLANT', 'YEAR'],
                    as_index=False)['VALUE'].sum()
    df = df.replace(color_dict)
    df = df.sort_values(by=['YEAR',
                            'NODE',
                            'POWERPLANT'])
    df['VALUE'] = df['VALUE'].round(2)
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
    df_capacities = powerplant_summary(df_capacities)

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
                                             'ProductionByTechnologyAnnual.csv'
                                             )
                                )
    df_generation = powerplant_summary(df_generation)

    return df_generation.to_csv(os.path.join(scenario_result_summaries_dir,
                                             'Generation.csv'
                                             ),
                                index=None
                                )


if __name__ == '__main__':
    main()
