import requests
import os
import yaml
import pandas as pd
from OPG_configuration import ConfigFile, ConfigPaths

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

    df_min_cap_inv = pd.read_csv(os.path.join(output_data_dir, 'TotalAnnualMinCapacityInvestment.csv'))
    df_min_cap_inv = df_min_cap_inv.append(tech_capacity_df)
    
    for each_year in list(years):
        
        
    
    df_max_cap = pd.read_csv(os.path.join(output_data_dir, 'TotalAnnualMaxCapacity.csv'))
    
    return print(tech_capacity_df.head(),
                 '\n',
                 df_min_cap_inv.head(),
                 '\n',
                 df_max_cap.head())


if __name__ == '__main__':
    main()
    logging.info(f'Max capacity limits sucessfully set')
