import os

# REQUIRED 

configfile: 'config/config.yaml'
output_dir = config['outputDir']
scenario = config['scenario']

# OUTPUT FILES 

result_files = [
    'AccumulatedNewCapacity.csv',
    'AnnualEmissions.csv',
    'AnnualFixedOperatingCost.csv',
    'AnnualTechnologyEmission.csv',
    'AnnualTechnologyEmissionByMode.csv',
    'AnnualVariableOperatingCost.csv',
    'CapitalInvestment.csv',
    'Demand.csv',
    'DiscountedTechnologyEmissionsPenalty.csv',
    'NewCapacity.csv',
    'ProductionByTechnology.csv',
    'ProductionByTechnologyAnnual.csv',
    'RateOfActivity.csv',
    'RateOfProductionByTechnology.csv',
    'RateOfProductionByTechnologyByMode.csv',
    'RateOfUseByTechnology.csv',
    'RateOfUseByTechnologyByMode.csv',
    'TotalAnnualTechnologyActivityByMode.csv',
    'TotalCapacityAnnual.csv',
    'TotalTechnologyAnnualActivity.csv',
    'TotalTechnologyModelPeriodActivity.csv',
    'UseByTechnology.csv'
]

result_figures = [
    'TotalCapacityAnnual', 
    'GenerationAnnual',
]

# RULES

rule otoole_results:
    message:
        'Generating result csv files...'
    input:
        solution_file = Path(output_dir, scenario, f'{scenario}.sol'),
        pre_process_file = Path(output_dir, scenario, f'PreProcessed_{scenario}.txt')
    output:
        expand(Path(output_dir, scenario, 'results/{result_file}'), result_file = result_files),
    conda:
        '../envs/otoole.yaml'
    log:
        'workflow/logs/otoole_results.log'
    shell: 
        'otoole results cbc csv {input.solution_file} {output_dir}/{scenario}/results ' 
        '--input_datafile {input.pre_process_file} '
        '--input_datapackage {output_dir}/{scenario}/datapackage.json ' 
        '2> {log}'

rule visualisation:
    message:
        'Generating result figures...'
    input:
        expand(Path(output_dir, scenario, 'results/{result_file}'), result_file = result_files),
        'config/config.yaml',
    output:
        expand(Path(output_dir, scenario, 'figures/{result_figure}.html'), result_figure = result_figures),
    conda:
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/visualisation.log'
    shell: 
        'python workflow/scripts/osemosys_global/visualisation.py 2> {log}'