'''Functionality to interface with configuration files. '''

from pathlib import Path
import yaml

class ConfigFile:
    '''Class to hold yaml configuration file data
    
    Args: 
        config_file_name = yaml file name in the config/ folder 
        
    Example: 
        config = ConfigFile('settings')
        config.get('geographic_scope')
        -> ['IND','NPL']
    '''

    # non changing parameters 
    region_name = 'GLOBAL' 

    def __init__ (self, config_file_name): 
        self.file_path = Path(Path(__file__).resolve().parent, 
            '../../../config', f'{config_file_name}.yaml')

    def get(self, name):
        with open(self.file_path, encoding='utf-8') as yaml_file:
            parsed_yaml_file = yaml.load(yaml_file, Loader = yaml.FullLoader).get(name)
        return parsed_yaml_file

    def get_years(self):
        start_year = self.get('startYear')
        end_year = self.get('endYear')
        return list(range(start_year, end_year + 1))

class ConfigPaths:
    '''Class to hold relative paths from file called from. '''    

    # Hard coded file structure 
    input_dir_name = 'resources'
    output_dir_name = 'results'
    py_file_dir = Path(__file__).resolve().parent # folder of this module 
    
    def __init__(self):
        self.input_dir = Path(self.py_file_dir, '../../../', self.input_dir_name)
        self.input_data_dir = Path(self.input_dir, 'data')

        self.output_dir = Path(self.py_file_dir, '../../../', self.output_dir_name)
        self.output_data_dir = Path(self.output_dir, 'data')

        self.scenario_dir = Path(self.output_dir, self.get_scenario_name())
        self.scenario_data_dir = Path(self.scenario_dir, 'data')
        self.scenario_figs_dir = Path(self.scenario_dir, 'figures')
        self.scenario_results_dir = Path(self.scenario_dir, 'results')
        self.scenario_result_summaries_dir = Path(self.scenario_dir,
                                                  'result_summaries')

        self.otoole = Path(self.input_dir, 'otoole')

        self.custom_nodes_dir = Path(self.input_dir, 'data/custom_nodes')

    def get_scenario_name(self):
        config = ConfigFile('config')
        return config.get('scenario') 