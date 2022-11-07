#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import sys
from OPG_configuration import ConfigPaths
from otoole.utils import _read_file
from pathlib import Path

def main(user_config):
    config_paths = ConfigPaths()
    output_data_dir = config_paths.output_data_dir

    for param, param_data in user_config.items():
        csv = f'{param}.csv'
        if csv in os.listdir(output_data_dir):
            continue
        if param_data['type'] == 'param':
            columns = param_data['indices']
            columns.append('VALUE')
            df = pd.DataFrame(columns=columns)
            df.to_csv(Path(output_data_dir, csv), index=False)
        if param_data['type'] == 'set':
            df = pd.DataFrame(columns=['VALUE'])
            df.to_csv(Path(output_data_dir, csv), index=False)

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python OPG_file_check.py <user_config>")
    else:
        _, ending = os.path.splitext(sys.argv[1])
        with open(sys.argv[1], "r") as f:
            user_config = _read_file(f, ending)

        main(user_config)