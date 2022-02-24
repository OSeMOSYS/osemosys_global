import os

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
    'variable_costs',
    'emissions'
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

# DATA PROCESSING RULES 

rule powerplant:
    message:
        'Generating powerplant data...'
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
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/powerplant.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_powerplant_data.py 2> {log}'

rule timeslice:
    message:
        'Generating timeslice data...'
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
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/timeslice.log'    
    shell:
        'python workflow/scripts/osemosys_global/OPG_TS_data.py 2> {log}'

rule variable_costs:
    message:
        'Generating variable cost data...'
    input:
        Path(input_dir, 'data/CMO-April-2020-forecasts.xlsx'),
        Path(output_dir, 'data/TECHNOLOGY.csv'),
        'config/config.yaml'
    output:
        expand(Path(output_dir, 'data/{output_file}'), output_file=variable_cost_files),
        touch('workflow/rules/flags/variable_costs.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/variable_costs.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_variablecosts.py 2> {log}'

rule demand_projections:
    message:
        'Generating demand data...'
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
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/demand_projections.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_demand_projection.py 2> {log}'

rule emissions:
    message:
        'Generating emission data...'
    input:
        Path(input_dir, 'data/emission_factors.csv'),
        Path(output_dir, 'data/InputActivityRatio.csv'),
        'config/config.yaml'
    output: 
        expand(Path(output_dir, 'data/{output_file}'), output_file = emission_files),
        touch('workflow/rules/flags/emissions.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/emissions.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_emissions.py 2> {log}'

