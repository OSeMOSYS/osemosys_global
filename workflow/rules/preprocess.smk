import os

# REQUIRED 

configfile: 'config/config.yaml'

# MODEL AND SCENARIO OUTPUT FILES 

osemosys_files = os.listdir('resources/simplicity/data')
osemosys_files.remove('default_values.csv') #taken from input_dir

demand_figures = [
    'South America',
    'Oceania',
    'North America',
    'Europe',
    'Asia',
    'Africa'
]

# Can instead replace with rule order
flag_files = [
    'demand_projections',
    'timeslice',
    'variable_costs',
    'emissions',
    'max_capacity'
]

# SCRIPT OUTPUT FILES 

power_plant_files = [
    'CapitalCost.csv',
    'FixedCost.csv',
    'CapacityToActivityUnit.csv',
    'OperationalLife.csv',
    'TotalAnnualMaxCapacityInvestment.csv',
    'TotalTechnologyModelPeriodActivityUpperLimit.csv',
    'FUEL.csv',
    'InputActivityRatio.csv',
    'OutputActivityRatio.csv',
    'MODE_OF_OPERATION.csv',
    'REGION.csv',
    'ResidualCapacity.csv',
    'TECHNOLOGY.csv',
    'YEAR.csv'
    ]

timeslice_files = [
    'CapacityFactor.csv',
    'TIMESLICE.csv',
    'SpecifiedDemandProfile.csv',
    'YearSplit.csv'
    ]

variable_cost_files = [
    'VariableCost.csv'
    ]

demand_files = [
    'SpecifiedAnnualDemand.csv'
    ]

emission_files = [
    'EmissionActivityRatio.csv',
    'EmissionsPenalty.csv',
    'EMISSION.csv'
]

max_capacity_files = [
    'TotalAnnualMaxCapacity.csv'
]

# DATA PROCESSING RULES 

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
        config['crossborderTrade'],
        config['startYear'],
        config['endYear'],
        config['no_invest_technologies']
    output:
        expand('results/data/{output_file}', output_file = power_plant_files)
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/logs/powerplant.log'
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
        config['startYear'],
        config['endYear'],
        config['daytype'],
        config['dayparts'],
        config['seasons'],
    output:
        csv_files = expand('results/data/{output_file}', output_file=timeslice_files),
        flag_file = touch('workflow/rules/flags/timeslice.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/logs/timeslice.log'    
    shell:
        'python workflow/scripts/osemosys_global/OPG_TS_data.py 2> {log}'

rule variable_costs:
    message:
        'Generating variable cost data...'
    input:
        'resources/data/CMO-April-2020-forecasts.xlsx',
        'results/data/TECHNOLOGY.csv',
    params:
        config['startYear'],
        config['endYear'],
    output:
        csv_files = expand('results/data/{output_file}', output_file=variable_cost_files),
        flag_file = touch('workflow/rules/flags/variable_costs.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/logs/variable_costs.log'
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
        config['startYear'],
        config['endYear'],
    output:
        csv_files = expand('results/data/{output_file}', output_file = demand_files),
        figures = expand('results/figs/Demand projection {demand_figure}.jpg', demand_figure = demand_figures),
        flag_file = touch('workflow/rules/flags/demand_projections.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/logs/demand_projections.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_demand_projection.py 2> {log}'

rule emissions:
    message:
        'Generating emission data...'
    input:
        'resources/data/emission_factors.csv',
        'results/data/InputActivityRatio.csv',
    params:
        config['startYear'],
        config['endYear'],
        config['emission_penalty'],
    output: 
        csv_files = expand('results/data/{output_file}', output_file = emission_files),
        flag_file = touch('workflow/rules/flags/emissions.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/logs/emissions.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_emissions.py 2> {log}'

rule max_capacity:
    message: 
        'Generating capacity limits...'
    input:
        'resources/data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx',
        'results/data/ResidualCapacity.csv'
    params:
        config['startYear'],
        config['endYear'],
    output:
        csv_files = expand('results/data/{output_file}', output_file = max_capacity_files),
        flag_file = touch('workflow/rules/flags/max_capacity.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/logs/max_capacity.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_max_capacity.py 2> {log}'