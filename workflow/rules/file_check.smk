import os
import shutil
import pandas as pd

# REQUIRED 

configfile: 'config/config.yaml'
input_dir = config['inputDir']
output_dir = config['outputDir']

# HELPER FILES

include: 'preprocess.smk'

# GET LIST OF MISSING FILES 

missing_files = os.listdir(Path(input_dir, 'simplicity/data'))
missing_files.remove('default_values.csv')
for csv in os.listdir(Path(output_dir, 'data')):
    if csv == 'default_values.csv':
        continue
    missing_files.remove(csv)

# RULE

rule FileCheck:
    input:
        expand('workflow/rules/flags/{flag_file}.done', flag_file = flag_files),
        expand(Path(input_dir, 'simplicity/data/{missing_file}'), missing_file = missing_files),
        Path(input_dir, 'data/default_values.csv')
    output: 
        expand(Path(output_dir, 'data/{missing_file}'), missing_file = missing_files),
        Path(output_dir, 'data/default_values.csv')
    log: 
        'workflow/logs/fileCheck.log'
    run: 
        shutil.copy(Path(input_dir, 'data/default_values.csv'),
                    Path(output_dir, 'data/default_values.csv'))
        for each_csv in os.listdir(Path(output_dir, 'data')):
            # Refresh the default values csv from resources/
            if each_csv == "default_values.csv":
                continue
            if each_csv not in os.listdir(Path(output_dir, 'data')): 
                csv_df_in = pd.read_csv(Path(input_dir, 'simplicity/data', each_csv))
                csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
                csv_df_out.to_csv(Path(output_dir,'data', each_csv), index = None)

'''
def get_missing_files():
    missing = os.listdir(Path(input_dir, 'simplicity/data'))
    for csv in os.listdir(Path(output_dir, 'data')):
        missing.remove(csv)
    #csv_paths = []
    #for csv in missing: 
    #    csv_paths.append(Path(input_dir, 'simplicity/data', csv))
    return missing

rule FileCheck:
    input:
        expand('resources/simplicity/data/{missing_file}', missing_file = get_missing_files),
        #expand('workflow/rules/flags/{flag_file}.done', flag_file = flag_files),
        #expand(Path(input_dir, 'simplicity/data/{osemosys_file}'), osemosys_file = osemosys_files),
        #Path(input_dir, 'data/default_values.csv'),
    output: 
        expand('results/data/{missing_file}', missing_file = get_missing_files),
        #expand(Path(output_dir,'data/{osemosys_file}'), osemosys_file = osemosys_files),
        #expand('results/data/{missing_file}.csv', missing_file = get_missing_files),
        #Path(output_dir, 'data/default_values.csv')
    log: 
        'workflow/logs/fileCheck.log'
    run: 
        #for each_csv in os.listdir(Path(input_dir, 'simplicity/data')):
        for each_csv in os.listdir(Path(input_dir, 'simplicity/data')): 
            if each_csv == "default_values.csv":
                    shutil.copy(Path(input_dir, 'data', each_csv),
                                Path(output_dir, 'data', each_csv))
            if each_csv not in os.listdir(Path(output_dir, 'data')): 
                csv_df_in = pd.read_csv(Path(input_dir, 'simplicity/data', each_csv))
                csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
                csv_df_out.to_csv(Path(output_dir,'data', each_csv), index = None)
'''
'''
filecheck.py searches through the existing set of created files and only 
creates missing ones. The script is written to create all missing files at once, 
want to only call the rule once. Therefore, we should not be using wildcards. 
If we want to use the expand() function, we need the list of already created
files at the beginning of the workflow, which is currently implemented and 
kinda ugly... This rule should probably be revisitied to improve its ledgibility.
'''

