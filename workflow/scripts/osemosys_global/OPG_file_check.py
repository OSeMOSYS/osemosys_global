#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import shutil
import yaml
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


#Read in information from YAML file
yaml_file = open(os.path.join(os.path.dirname(__file__), '../../..',
                              'config/config.yaml'))
parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

input_dir = os.path.join(os.path.dirname(__file__), '../../..',
    parsed_yaml_file.get('inputDir'))
input_data_dir = os.path.join(input_dir, 'data')
output_dir = os.path.join(os.path.dirname(__file__), '../../..', 
    parsed_yaml_file.get('outputDir'))
output_data_dir = os.path.join(output_dir, 'data')

simplicity_data = os.path.join(input_dir, 'simplicity/data')

for each_csv in os.listdir(simplicity_data):
    #copy default values csv with data
    if each_csv == "default_values.csv":
            shutil.copy(os.path.join(input_data_dir, each_csv),
                        os.path.join(output_data_dir, each_csv))
    if each_csv not in os.listdir(output_data_dir):
        csv_df_in = pd.read_csv(os.path.join(simplicity_data, each_csv))
        csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
        csv_df_out.to_csv(os.path.join(output_data_dir, each_csv),
                             index = None)
    

logging.info('File Check Completed')
