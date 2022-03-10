#!/usr/bin/env python
# coding: utf-8

# # Filter osemosys_global datapackaged based on user-defined geographic scope

import pandas as pd
import os
import yaml
import shutil
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

_PY_DIR = os.path.dirname(__file__)
yaml_file = open(os.path.join(_PY_DIR, '../../../config/config.yaml'))
parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

scenario_name = parsed_yaml_file.get('scenario')
input_dir = os.path.join(_PY_DIR, '../../..', parsed_yaml_file.get('outputDir'), 'data')
output_dir = os.path.join(_PY_DIR, '../../..', parsed_yaml_file.get('outputDir'),  scenario_name, 'data')

geographic_scope = parsed_yaml_file.get('geographic_scope')
if not geographic_scope: # Check for empty list (ie. World run)
    geographic_scope = []
geographic_scope.append('INT') # 'INT' for international fuels added by default
international_fuels = ['COA', 'COG', 'GAS', 'OIL', 'PET', 'OTH', 'URN']

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for each_csv in (os.listdir(input_dir)):
    df = pd.read_csv(os.path.join(input_dir,each_csv))

    if not df.empty:
        # Do not filter if only element is international fuels
        if geographic_scope[0] != 'INT': 
            if 'TECHNOLOGY' in df.columns:
                df = df.loc[df['TECHNOLOGY'].str[3:6].isin(geographic_scope) | 
                            df['TECHNOLOGY'].str[6:9].isin(geographic_scope) | 
                            df['TECHNOLOGY'].str[8:11].isin(geographic_scope)]

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
        
    df.to_csv(os.path.join(os.path.join(output_dir,each_csv)), index = None)

# copy datapackage over for otoole convert
shutil.copyfile(os.path.join(_PY_DIR, '../../..',
                             parsed_yaml_file.get('inputDir'), 
                             'simplicity/datapackage.json'),
                os.path.join(_PY_DIR, '../../..',
                             parsed_yaml_file.get('outputDir'),
                             scenario_name,
                             'datapackage.json')
)

logging.info('Geographic Filter Applied')