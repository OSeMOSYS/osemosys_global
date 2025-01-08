import subprocess
import os
import shutil

config_dir = 'config_consecutive/3hourly2seasons'
data_dir = 'results/data'

'''Set to 'list' to only run scenarios in 'scenario_list' or set to 'folder' to run all
scenario config files as defined in 'config_dir'. Can be used to run a scenario sample.'''
run_type = 'list'

scenario_list = [
  'ASEAN',
  'MYSPESGPXX',
  'MYSPETHASO',
  'IDNSMMYSPE',
  'MYSSHPHLLU',
  'BRNXXMYSSK',
  'LAOXXTHANO',
  'LAOXXVNMNO',
  'MMRXXTHANO',
  'KHMXXLAOXX',
  'KHMXXTHACE',
  'IDNKAMYSSH',
  'IDNSMSGPXX',
  'LAOXXMMRXX',
  'SGPXXVNMSO',
  'KHMXXSGPXX',
  'MYSPEMYSSK',
  'MYSSHMYSSK',
  'IDNJWIDNKA',
  'IDNJWIDNSM',
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