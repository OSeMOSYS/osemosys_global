import pandas as pd
import itertools
import os
import sys
import yaml
from OPG_configuration import ConfigFile, ConfigPaths
pd.set_option('mode.chained_assignment', None)

def main():
    config_paths = ConfigPaths()
    config = ConfigFile('config')

    scenario_result_summaries_dir = config_paths.scenario_result_summaries_dir
    
    with open(scenario_result_summaries_dir)
    
if __name__ == '__main__':
    main()