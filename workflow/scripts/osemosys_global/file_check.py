#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import sys
from pathlib import Path
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def main(data: str, otoole_data: str):
    """Checks for missing otoole files"""
    DATA = Path(data)
    OTOOLE = Path(otoole_data)
    
    for each_csv in os.listdir(OTOOLE):
        if each_csv not in os.listdir(DATA): 
            csv_df_in = pd.read_csv(Path(OTOOLE, each_csv))
            csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
            csv_df_out.to_csv(Path(DATA, each_csv), index = None)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python file_check.py.py <input/csv/dir/> </otoole/csv/dir/> ")
        logging.info('File Check Failed')
    else:
        main(
            sys.argv[1], 
            sys.argv[2], 
        )
        logging.info('File Check Completed')
