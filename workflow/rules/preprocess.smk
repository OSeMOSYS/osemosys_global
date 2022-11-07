import os
import shutil
configfile: 'config/config.yaml'
include: 'data_files.smk'

rule make_data_dir:
    output: directory('results/data')
    shell: 'mkdir -p {output}'

rule powerplant:
    message:
        'Generating powerplant data...'
    input:
        'resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx',
        'resources/data/weo_2018_powerplant_costs.csv',
        'resources/data/operational_life.csv',
        'resources/data/naming_convention_tech.csv',
        'resources/data/Costs Line expansion.xlsx',
        'resources/data/weo_region_mapping.csv',
    params: 
        trade = config['crossborderTrade'],
        start_year = config['startYear'],
        end_year = config['endYear'],
        invest_techs = config['no_invest_technologies']
    output:
        csv_files = expand('results/data/{output_file}', output_file = POWER_PLANT_FILES)
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/data/logs/powerplant.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_powerplant_data.py 2> {log}'

rule timeslice:
    message:
        'Generating timeslice data...'
    input:
        'resources/data/All_Demand_UTC_2015.csv',
        'resources/data/CSP 2015.csv',
        'resources/data/SolarPV 2015.csv',
        'resources/data/Hydro_Monthly_Profiles (15 year average).csv',
        'resources/data/Won 2015.csv',
        'resources/data/Woff 2015.csv',
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        daytype = config['daytype'],
        daypart = config['dayparts'],
        seasons = config['seasons'],
    output:
        csv_files = expand('results/data/{output_file}', output_file=TIMESLICE_FILES),
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/data/logs/timeslice.log'    
    shell:
        'python workflow/scripts/osemosys_global/OPG_TS_data.py 2> {log}'

rule variable_costs:
    message:
        'Generating variable cost data...'
    input:
        'resources/data/CMO-April-2020-forecasts.xlsx',
        'results/data/TECHNOLOGY.csv',
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
    output:
        csv_files = expand('results/data/{output_file}', output_file=VARIABLE_COST_FILES),
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/data/logs/variable_costs.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_variablecosts.py 2> {log}'

rule demand_projections:
    message:
        'Generating demand data...'
    input:
        'resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx',
        'resources/data/iamc_db_GDPppp_Countries.xlsx',
        'resources/data/iamc_db_POP_Countries.xlsx',
        'resources/data/iamc_db_URB_Countries.xlsx',
        'resources/data/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx',
        'resources/data/T&D Losses.xlsx',
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
    output:
        csv_files = expand('results/data/{output_file}', output_file = DEMAND_FILES),
        figures = expand('results/data/../figs/Demand projection {demand_figure}.jpg', demand_figure = DEMAND_FIGURES),
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/data/logs/demand_projections.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_demand_projection.py 2> {log}'

rule emissions:
    message:
        'Generating emission data...'
    input:
        'resources/data/emission_factors.csv',
        'results/data/InputActivityRatio.csv',
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        emission = config['emission_penalty']
    output: 
        csv_files = expand('results/data/{output_file}', output_file = EMISSION_FILES),
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/data/../logs/emissions.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_emissions.py 2> {log}'

rule max_capacity:
    message: 
        'Generating capacity limits...'
    input:
        'resources/data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx',
        'results/data/ResidualCapacity.csv'
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
    output:
        csv_files = expand('results/data/{output_file}', output_file = MAX_CAPACITY_FILES),
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/data/logs/max_capacity.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_max_capacity.py 2> {log}'

rule copy_default_values:
    message:
        'Copying user defined default values...'
    input:
        'resources/data/default_values.csv'
    output:
        'results/data/default_values.csv'
    run:
        shutil.copyfile(input[0], output[0])

rule file_check:
    message:
        'Generating missing files...'
    input:
        csv_files = expand('results/data/{output_file}', output_file = GENERATED_FILES),
    params:
        config = 'resources/config.yaml'
    output: 
        expand('results/data/{check_file}', check_file = CHECK_FILES),
    conda:
        '../envs/data_processing.yaml'
    log: 
        'results/data/logs/file_check.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_file_check.py {params.config} 2> {log}'


