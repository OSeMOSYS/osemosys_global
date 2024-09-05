# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:14:43 2024

@author: maart
"""

import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
from pathlib import Path
from configuration import ConfigPaths
import os
import requests

# CONFIGURATION PARAMETERS
config_paths = ConfigPaths()
input_data_dir = config_paths.input_data_dir

external_files = {
    'PLEXOS_World_2015_Gold_V1.1.xlsx' : 
    'https://dataverse.harvard.edu/api/access/datafile/4008393?format=original&gbrecs=true',
    
    'All_Demand_UTC_2015.csv' :
    'https://dataverse.harvard.edu/api/access/datafile/3985039?format=original&gbrecs=true',
    
    'PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx' :
    'https://dataverse.harvard.edu/api/access/datafile/6040815'
    
                  }

if __name__ == "__main__":
    for file, url in external_files.items():
        path = os.path.join(input_data_dir, file)

        if not Path(path).exists():
            logging.info(f'Downloading {file}')

            data = requests.get(url , path)
            
            with open(path, 'wb') as f:
                f.write(data.content)