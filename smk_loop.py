import subprocess
import os
import shutil

config_dir = 'config_consecutive'
data_dir = 'results/data'

'''Set to 'list' to only run scenarios in 'scenario_list' or set to 'folder' to run all
scenario config files as defined in 'config_dir'. Can be used to run a scenario sample.'''
run_type = 'folder'

scenario_list = [
    'lowseasons',
    'base',
    'trnlimit',
    'reservemargin',
    'varcost',
    #'trn20',
    #'trn21',
    #'trn22',
    #'trn23',
    #'trn24',
    #'trn25',
    #'trn26',
    #'trn27',
    #'trn28',
    #'trn29',
    #'trn30',
    #'trn31',
    #'trn32',
    #'trn33',
    #'trn34',
    ]

if run_type == 'folder':
    for scenario in os.listdir(config_dir):
        subprocess.run("snakemake --cores 12 --configfile "f'{config_dir}/{scenario}'"", 
                       shell = True)
        
        shutil.rmtree(data_dir)
        
if run_type == 'list':
    for scenario in scenario_list:
        subprocess.run("snakemake --cores 12 --configfile "f'{config_dir}/{scenario}.yaml'"", 
                       shell = True)
        
        shutil.rmtree(data_dir)