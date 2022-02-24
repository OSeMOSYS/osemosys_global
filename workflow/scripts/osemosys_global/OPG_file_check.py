#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import shutil
import yaml
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# File Location

_PY_DIR = os.path.dirname(__file__)

#Read in information from YAML file

yaml_file = open(os.path.join(_PY_DIR, '../../..', 'config/config.yaml'))
parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

# Standard Paths

input_dir = os.path.join(_PY_DIR, '../../..', parsed_yaml_file.get('inputDir'))
output_dir = os.path.join(_PY_DIR, '../../..', parsed_yaml_file.get('outputDir'))
simplicity_data = os.path.join(input_dir, 'simplicity/data')

# File Check logic

shutil.copy(os.path.join(input_dir, 'data/default_values.csv'),
            os.path.join(output_dir, 'data/default_values.csv'))
for each_csv in os.listdir(simplicity_data):
    # Default values csv already copied from resources/
    if each_csv == "default_values.csv":
        continue
    if each_csv not in os.listdir(os.path.join(output_dir, 'data')): 
        csv_df_in = pd.read_csv(os.path.join(simplicity_data, each_csv))
        csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
        csv_df_out.to_csv(os.path.join(output_dir,'data', each_csv), index = None)
    
logging.info('File Check Completed')
