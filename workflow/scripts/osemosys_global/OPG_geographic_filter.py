#!/usr/bin/env python
# coding: utf-8

# # Filter osemosys_global datapackaged based on user-defined geographic scope

import pandas as pd
import os
import yaml
import shutil
from pathlib import Path
from OPG_configuration import ConfigFile, ConfigPaths
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# CONFIGURATION PARAMETERS

config_paths = ConfigPaths()
config = ConfigFile('config')  

scenario_name = config.get('scenario')
geographic_scope = config.get('geographic_scope')
remove_nodes = config.get('nodes_to_remove')

output_data_dir = config_paths.output_data_dir
scenario_dir = config_paths.scenario_dir
scenario_data_dir = config_paths.scenario_data_dir

# FILTERING 

if not geographic_scope: # Check for empty list (ie. World run)
    geographic_scope = []
geographic_scope.append('INT') # 'INT' for international fuels added by default
international_fuels = ['COA', 'COG', 'GAS', 'OIL', 'PET', 'OTH', 'URN']

if not os.path.exists(scenario_data_dir):
    os.makedirs(scenario_data_dir)

for each_csv in Path(output_data_dir).glob('*.csv'):
    df = pd.read_csv(os.path.join(output_data_dir, each_csv))

    if not df.empty:
        # Do not filter if only element is international fuels
        if geographic_scope[0] != 'INT': 
            if 'TECHNOLOGY' in df.columns:
                df = df.loc[df['TECHNOLOGY'].str[3:6].isin(geographic_scope) | 
                            df['TECHNOLOGY'].str[6:9].isin(geographic_scope) | 
                            df['TECHNOLOGY'].str[8:11].isin(geographic_scope)]

                # Filter out all international TRN techs 
                df = df.loc[~(
                    df['TECHNOLOGY'].str.startswith('TRN') &
                    (~(df['TECHNOLOGY'].str[3:6].isin(geographic_scope)) |
                    ~(df['TECHNOLOGY'].str[8:11].isin(geographic_scope)))
                    )]

                if remove_nodes:
                    df = df.loc[~(df['TECHNOLOGY'].str[3:8].isin(remove_nodes) | 
                                df['TECHNOLOGY'].str[6:11].isin(remove_nodes) | 
                                df['TECHNOLOGY'].str[8:13].isin(remove_nodes))]

            if 'FUEL' in df.columns:
                df = df.loc[df['FUEL'].str[3:6].isin(geographic_scope) | 
                            df['FUEL'].str[6:9].isin(geographic_scope) |
                            df['FUEL'].isin(international_fuels)]
                
                if remove_nodes:    
                    df = df.loc[~(df['FUEL'].str[3:8].isin(remove_nodes) | 
                                df['FUEL'].str[6:11].isin(remove_nodes))]

            if str(each_csv).split('/')[-1] == 'FUEL.csv':
                df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                            df['VALUE'].str[6:9].isin(geographic_scope) |
                            df['VALUE'].isin(international_fuels)]
                
                if remove_nodes:
                    df = df.loc[~(df['VALUE'].str[3:8].isin(remove_nodes) | 
                                df['VALUE'].str[6:11].isin(remove_nodes))]

            if str(each_csv).split('/')[-1] == 'TECHNOLOGY.csv':
                df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                            df['VALUE'].str[6:9].isin(geographic_scope) | 
                            df['VALUE'].str[8:11].isin(geographic_scope)]
                df = df.loc[~(
                    df['VALUE'].str.startswith('TRN') &
                    (~(df['VALUE'].str[3:6].isin(geographic_scope)) |
                    ~(df['VALUE'].str[8:11].isin(geographic_scope)))
                    )]

                if remove_nodes:
                    df = df.loc[~(df['VALUE'].str[3:8].isin(remove_nodes) | 
                                df['VALUE'].str[6:11].isin(remove_nodes) | 
                                df['VALUE'].str[8:13].isin(remove_nodes))]
        
    df.to_csv(os.path.join(os.path.join(scenario_data_dir, each_csv.name)), index = None)


logging.info('Geographic Filter Applied')
