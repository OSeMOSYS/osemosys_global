#!/usr/bin/env python
# coding: utf-8

# # Filter osemosys_global datapackaged based on user-defined geographic scope

import pandas as pd
import os
import yaml
import shutil
from OPG_configuration import ConfigFile, ConfigPaths
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# CONFIGURATION PARAMETERS

config_paths = ConfigPaths()
config = ConfigFile('config')  

scenario_name = config.get('scenario')
geographic_scope = config.get('geographic_scope')

output_data_dir = config_paths.output_data_dir
scenario_dir = config_paths.scenario_dir
scenario_data_dir = config_paths.scenario_data_dir
simplicity_dir = config_paths.simplicity

# FILTERING 

if not geographic_scope: # Check for empty list (ie. World run)
    geographic_scope = []
geographic_scope.append('INT') # 'INT' for international fuels added by default
international_fuels = ['COA', 'COG', 'GAS', 'OIL', 'PET', 'OTH', 'URN']

if not os.path.exists(scenario_data_dir):
    os.makedirs(scenario_data_dir)

for each_csv in (os.listdir(output_data_dir)):
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

            if 'FUEL' in df.columns:
                df = df.loc[df['FUEL'].str[3:6].isin(geographic_scope) | 
                            df['FUEL'].str[6:9].isin(geographic_scope) |
                            df['FUEL'].isin(international_fuels)]

            if each_csv == 'FUEL.csv':
                df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                            df['VALUE'].str[6:9].isin(geographic_scope) |
                            df['VALUE'].isin(international_fuels)]

            if each_csv == 'TECHNOLOGY.csv':
                df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                            df['VALUE'].str[6:9].isin(geographic_scope) | 
                            df['VALUE'].str[8:11].isin(geographic_scope)]

                df = df.loc[~(
                    df['VALUE'].str.startswith('TRN') &
                    (~(df['VALUE'].str[3:6].isin(geographic_scope)) |
                    ~(df['VALUE'].str[8:11].isin(geographic_scope)))
                    )]
        
    df.to_csv(os.path.join(os.path.join(scenario_data_dir, each_csv)), index = None)

# copy datapackage over for otoole convert
shutil.copyfile(os.path.join(simplicity_dir, 'datapackage.json'),
                os.path.join(scenario_dir, 'datapackage.json'))

logging.info('Geographic Filter Applied')
