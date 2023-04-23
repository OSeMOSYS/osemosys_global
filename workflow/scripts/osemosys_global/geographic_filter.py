#!/usr/bin/env python
# coding: utf-8

# # Filter osemosys_global datapackaged based on user-defined geographic scope

import pandas as pd
import os
import sys
from pathlib import Path
from typing import List
from osemosys_global.configuration import ConfigFile, ConfigPaths
from osemosys_global.utils import get_config_data
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def main(input_data: str, output_data: str, geographic_scope: List[str], nodes_to_remove: List[str]):
    """Filters for scenario geography 
    
    Args:
        input_data: str
            Path to resources/data folder 
        output_dir: str
            Path to results/data folder
        geographic_scope: List[str]
            Start year of demand 
        nodes_to_remove: List[str]
            Nodes to remove from the model
    """
    
    # set path variables
    INPUT_DATA_DIR = Path(input_data)
    OUTPUT_DATA_DIR = Path(output_data)
    
    if not geographic_scope: # Check for empty list (ie. World run)
        geographic_scope = []
    geographic_scope.append('INT') # 'INT' for international fuels added by default
    international_fuels = ['COA', 'COG', 'GAS', 'OIL', 'PET', 'OTH', 'URN']

    if not os.path.exists(OUTPUT_DATA_DIR):
        os.makedirs(OUTPUT_DATA_DIR)

    for each_csv in Path(INPUT_DATA_DIR).glob('*.csv'):
        df = pd.read_csv(os.path.join(INPUT_DATA_DIR, each_csv))

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

                    if nodes_to_remove:
                        df = df.loc[~(df['TECHNOLOGY'].str[3:8].isin(nodes_to_remove) | 
                                    df['TECHNOLOGY'].str[6:11].isin(nodes_to_remove) | 
                                    df['TECHNOLOGY'].str[8:13].isin(nodes_to_remove))]

                if 'FUEL' in df.columns:
                    df = df.loc[df['FUEL'].str[3:6].isin(geographic_scope) | 
                                df['FUEL'].str[6:9].isin(geographic_scope) |
                                df['FUEL'].isin(international_fuels)]
                    
                    if nodes_to_remove:    
                        df = df.loc[~(df['FUEL'].str[3:8].isin(nodes_to_remove) | 
                                    df['FUEL'].str[6:11].isin(nodes_to_remove))]

                if str(each_csv).split('/')[-1] == 'FUEL.csv':
                    df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                                df['VALUE'].str[6:9].isin(geographic_scope) |
                                df['VALUE'].isin(international_fuels)]
                    
                    if nodes_to_remove:
                        df = df.loc[~(df['VALUE'].str[3:8].isin(nodes_to_remove) | 
                                    df['VALUE'].str[6:11].isin(nodes_to_remove))]

                if str(each_csv).split('/')[-1] == 'TECHNOLOGY.csv':
                    df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                                df['VALUE'].str[6:9].isin(geographic_scope) | 
                                df['VALUE'].str[8:11].isin(geographic_scope)]
                    df = df.loc[~(
                        df['VALUE'].str.startswith('TRN') &
                        (~(df['VALUE'].str[3:6].isin(geographic_scope)) |
                        ~(df['VALUE'].str[8:11].isin(geographic_scope)))
                        )]

                    if nodes_to_remove:
                        df = df.loc[~(df['VALUE'].str[3:8].isin(nodes_to_remove) | 
                                    df['VALUE'].str[6:11].isin(nodes_to_remove) | 
                                    df['VALUE'].str[8:13].isin(nodes_to_remove))]
            
        df.to_csv(os.path.join(os.path.join(scenario_data_dir, each_csv.name)), index = None)
        df.to_csv(Path(OUTPUT_DIR, each_csv.name), index = None)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python geographic_filter.py <input/csv/dir> <output/csv/dir/> <config.yaml>")
        logging.info('Geographic Filter Completed')
    else:
        config_values = ["geographic_scope", "nodes_to_remove"]
        config_data = get_config_data(sys.argv[3], config_values)
        
        main(
            sys.argv[1], 
            sys.argv[2], 
            config_data["geographic_scope"],
            config_data["nodes_to_remove"],
        )
        logging.info('Geographic Filter Completed')
