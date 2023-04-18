#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import sys
from osemosys_global.configuration import ConfigPaths
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def main(output_data: str, otoole_data: str):
    """Checks for missing otoole files"""
    for each_csv in os.listdir(otoole_data):
        if each_csv not in os.listdir(output_data): 
            csv_df_in = pd.read_csv(os.path.join(otoole_data, each_csv))
            csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
            csv_df_out.to_csv(os.path.join(output_data, each_csv), index = None)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python demand_projection.py <results/data> <resources/otoole> ")
        logging.info('File Check Failed')
    else:
        main(
            sys.argv[1], 
            sys.argv[2], 
        )
        logging.info('File Check Completed')
