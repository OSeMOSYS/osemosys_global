#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import shutil
# from osemosys_global.configuration import ConfigPaths
from configuration import ConfigFile, ConfigPaths
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


# CONFIGURATION PARAMETERS
config_paths = ConfigPaths()
input_data_dir = config_paths.input_data_dir
output_data_dir = config_paths.output_data_dir
otoole_data = os.path.join(config_paths.otoole, "data")

# File Check logic

for each_csv in os.listdir(otoole_data):    
    if each_csv not in os.listdir(output_data_dir): 
        csv_df_in = pd.read_csv(os.path.join(otoole_data, each_csv))
        csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
        csv_df_out.to_csv(os.path.join(output_data_dir, each_csv), index = None)
    
logging.info('File Check Completed')
