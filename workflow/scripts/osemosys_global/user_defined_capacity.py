
import os
import pandas as pd
from osemosys_global.OPG_configuration import ConfigFile, ConfigPaths
from utils import apply_dtypes

# LOGGING
import logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def main():
    '''Creates capacity limits on renewable technologies.'''
    
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile('config')  

    input_dir = config_paths.input_dir
    output_data_dir = config_paths.output_data_dir
    scenario_data_dir = config_paths.scenario_data_dir
    region = config.get('region')
    years = range(config.get('startYear'), config.get('endYear') + 1)
    tech_capacity = config.get('user_defined_capacity')
    techCapacity = []
    
    for tech, tech_params in tech_capacity.items():
        techCapacity.append([tech, tech_params[0], tech_params[1]])
    tech_capacity_df = pd.DataFrame(techCapacity,
                                    columns=['TECHNOLOGY', 'VALUE', 'YEAR'])
    tech_capacity_df['REGION'] = region
    tech_capacity_df = tech_capacity_df[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]

    tech_set = pd.read_csv(os.path.join(scenario_data_dir, 'TECHNOLOGY.csv'))

    for each_tech in list(tech_capacity_df['TECHNOLOGY'].unique()):
        if each_tech not in list(tech_set['VALUE']):
            tech_capacity_df = tech_capacity_df.loc[~(tech_capacity_df['TECHNOLOGY'].isin([each_tech]))]
    df_min_cap_inv = pd.read_csv(os.path.join(scenario_data_dir, 'TotalAnnualMinCapacityInvestment.csv'))
    df_min_cap_inv = df_min_cap_inv.append(tech_capacity_df)
    df_min_cap_inv.drop_duplicates(inplace=True)

    df_max_cap_inv = pd.read_csv(os.path.join(scenario_data_dir, 'TotalAnnualMaxCapacityInvestment.csv'))

    max_cap_techs = []
    for index, row in tech_capacity_df.iterrows():
        for each_year in years:
            if row['YEAR'] == each_year:
                value = row['VALUE']
            else:
                value = 0
            max_cap_techs.append([row['REGION'],
                                  row['TECHNOLOGY'],
                                  each_year,
                                  value])
    max_cap_techs_df = pd.DataFrame(max_cap_techs,
                                    columns=['REGION',
                                             'TECHNOLOGY',
                                             'YEAR',
                                             'VALUE'])
    df_max_cap_inv = df_max_cap_inv.append(max_cap_techs_df)
    df_max_cap_inv.drop_duplicates(inplace=True)

    df_max_cap_inv = apply_dtypes(df_max_cap_inv, "TotalAnnualMaxCapacityInvestment")
    df_min_cap_inv = apply_dtypes(df_min_cap_inv, "TotalAnnualMinCapacityInvestment")

    df_max_cap_inv.to_csv(os.path.join(
        scenario_data_dir, "TotalAnnualMaxCapacityInvestment.csv"), index=None)
    df_min_cap_inv.to_csv(os.path.join(
        scenario_data_dir, "TotalAnnualMinCapacityInvestment.csv"), index=None)

if __name__ == '__main__':
    main()
    logging.info(f'User-defined capacities sucessfully set')
