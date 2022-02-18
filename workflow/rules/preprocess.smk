import os
import shutil
import pandas as pd

# REQUIRED 

configfile: 'config/config.yaml'
input_dir = config['inputDir']
output_dir = config['outputDir']

# MODEL AND SCENARIO OUTPUT FILES 

osemosys_files = os.listdir(Path(input_dir, 'simplicity/data'))
osemosys_files.remove('default_values.csv') #taken from input_dir

demand_figures = [
    'South America',
    'Oceania',
    'North America',
    'Europe',
    'Asia',
    'Africa'
]

flag_files = [
    'demand_projections',
    'timeslice',
    'variable_costs'
]

# OUTPUT FILES FROM SCRIPTS 

# Files to get created with empty data
default_files = os.listdir(Path(input_dir, 'simplicity/data')) 

power_plant_files = [
    'CapitalCost.csv',
    'FixedCost.csv',
    'CapacityToActivityUnit.csv',
    'OperationalLife.csv',
    'TotalAnnualMaxCapacityInvestment.csv',
    'EMISSION.csv',
    'FUEL.csv',
    'InputActivityRatio.csv',
    'MODE_OF_OPERATION.csv',
    'OutputActivityRatio.csv',
    'REGION.csv',
    'ResidualCapacity.csv',
    'TECHNOLOGY.csv',
    'YEAR.csv'
    ]
[default_files.remove(csv) for csv in power_plant_files]

timeslice_files = [
    'CapacityFactor.csv',
    'TIMESLICE.csv',
    'SpecifiedDemandProfile.csv',
    'YearSplit.csv'
    ]
[default_files.remove(csv) for csv in timeslice_files]

variable_cost_files = [
    'VariableCost.csv'
    ]
[default_files.remove(csv) for csv in variable_cost_files]

demand_files = [
    'SpecifiedAnnualDemand.csv'
    ]
[default_files.remove(csv) for csv in demand_files]

# DATA PROCESSING RULES 

rule PowerPlant:
    input:
        Path(input_dir, 'data/PLEXOS_World_2015_Gold_V1.1.xlsx'),
        Path(input_dir, 'data/weo_2018_powerplant_costs.csv'),
        Path(input_dir, 'data/operational_life.csv'),
        Path(input_dir, 'data/naming_convention_tech.csv'),
        Path(input_dir, 'data/Costs Line expansion.xlsx'),
        Path(input_dir, 'data/weo_region_mapping.csv'),
        'config/config.yaml'
    output:
        expand(Path(output_dir, 'data/{output_file}'), output_file = power_plant_files)
    conda:
        'envs/powerPlant.yaml'
    log:
        'workflow/logs/PowerPlant.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_powerplant_data.py 2> {log}'

rule TimeSlice:
    input:
        Path(input_dir, 'data/All_Demand_UTC_2015.csv'),
        Path(input_dir, 'data/CSP 2015.csv'),
        Path(input_dir, 'data/SolarPV 2015.csv'),
        Path(input_dir, 'data/Hydro_Monthly_Profiles (15 year average).csv'),
        Path(input_dir, 'data/Won 2015.csv'),
        Path(input_dir, 'data/Woff 2015.csv'),
        'config/config.yaml'
    output:
        expand(Path(output_dir, 'data/{output_file}'), output_file=timeslice_files),
        touch('workflow/rules/flags/timeslice.done')
    conda:
        'envs/timeSlice.yaml'
    log:
        'workflow/logs/timeSlice.log'    
    shell:
        'python workflow/scripts/osemosys_global/OPG_TS_data.py 2> {log}'

rule VariableCosts:
    input:
        Path(input_dir, 'data/CMO-April-2020-forecasts.xlsx'),
        Path(output_dir, 'data/TECHNOLOGY.csv'),
        'config/config.yaml'
    output:
        expand(Path(output_dir, 'data/{output_file}'), output_file=variable_cost_files),
        touch('workflow/rules/flags/variable_costs.done')
    conda:
        'envs/variableCosts.yaml'
    log:
        'workflow/logs/variableCosts.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_variablecosts.py 2> {log}'

rule DemandProjections:
    input:
        Path(input_dir, 'data/PLEXOS_World_2015_Gold_V1.1.xlsx'),
        Path(input_dir, 'data/iamc_db_GDPppp_Countries.xlsx'),
        Path(input_dir, 'data/iamc_db_POP_Countries.xlsx'),
        Path(input_dir, 'data/iamc_db_URB_Countries.xlsx'),
        Path(input_dir, 'data/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx'),
        Path(input_dir, 'data/T&D Losses.xlsx'),
        'config/config.yaml'
    output:
        expand(Path(output_dir, 'data/{output_file}'), output_file = demand_files),
        expand(Path(output_dir, 'figs/Demand projection {demand_figure}.jpg'), demand_figure = demand_figures),
        touch('workflow/rules/flags/demand_projections.done')
    conda:
        'envs/demand.yaml'
    log:
        'workflow/logs/demandProjections.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_demand_projection.py 2> {log}'

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


filecheck.py searches through the existing set of created files and only 
creates missing ones. The script is written to create all missing files at once, 
want to only call the rule once. Therefore, we should not be using wildcards. 
If we want to use the expand() function, we need the list of already created
files at the beginning of the workflow, which is currently implemented and 
kinda ugly... This rule should probably be revisitied to improve its ledgibility.

rule FileCheck:
    input:
        expand('workflow/rules/flags/{flag_file}.done', flag_file = flag_files),
        expand(Path(input_dir, 'simplicity/data/{default_file}'), default_file = default_files),
        Path(input_dir, 'data/default_values.csv')
    output: 
        expand(Path(output_dir, 'data/{default_file}'), default_file = default_files),
    log: 
        'workflow/logs/fileCheck.log'
    run: 
        for each_csv in default_files:
            print(each_csv)
            if each_csv == "default_values.csv":
                    shutil.copy(Path(input_dir, 'data', each_csv),
                                Path(output_dir, 'data', each_csv))
            #if each_csv not in os.listdir(Path(output_dir, 'data')): 
            else:
                csv_df_in = pd.read_csv(Path(input_dir, 'simplicity/data', each_csv))
                csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
                csv_df_out.to_csv(Path(output_dir,'data', each_csv), index = None)



''' 





