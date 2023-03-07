import os

# REQUIRED 

configfile: 'config/config.yaml'

# OUTPUT FILES 

result_files = [
    'AccumulatedNewCapacity.csv',
    'AnnualEmissions.csv',
    'AnnualFixedOperatingCost.csv',
    'AnnualTechnologyEmission.csv',
    'AnnualTechnologyEmissionByMode.csv',
    'AnnualVariableOperatingCost.csv',
    # 'CapitalInvestment.csv',
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

result_summaries = [
    'Metrics'
]

# imput functions

def solver_file_type(wildcards):
    if config['solver'] == 'cplex':
        return 'results/{scenario}/{scenario}_sort.sol'
    else: 
        return 'results/{scenario}/{scenario}.sol'

# rules

rule otoole_results:
    message:
        'Generating result csv files...'
    input:
        solution_file = solver_file_type,
        pre_process_file = 'results/{scenario}/{scenario}.txt',
        otoole_config = 'results/{scenario}/otoole.yaml',
    output:
        expand('results/{{scenario}}/results/{result_file}', result_file = result_files),
    log:
        log = 'results/{scenario}/logs/otoole_results.log'
    shell: 
        '''
        otoole results {config[solver]} csv \
        {input.solution_file} results/{wildcards.scenario}/results \
        {input.otoole_config} --input_datafile {input.pre_process_file} \
        2> {log} 
        '''

rule visualisation:
    message:
        'Generating result figures...'
    input:
        csv_files = expand('results/{{scenario}}/results/{result_file}', result_file = result_files),
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        dayparts = config['dayparts'],
        seasons = config['seasons'],
        by_country = config['results_by_country'],
        geographic_scope = config['geographic_scope'],
    output:
        expand('results/{{scenario}}/figures/{result_figure}.html', result_figure = result_figures)
    log:
        log = 'results/{scenario}/logs/visualisation.log'
    shell: 
        'python workflow/scripts/osemosys_global/visualise.py 2> {log}'

rule summarise_results:
    message:
        'Generating summary of results...'
    input:
        csv_files = expand('results/{{scenario}}/results/{result_file}', result_file = result_files),
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        dayparts = config['dayparts'],
        seasons = config['seasons'],
    output:
        expand('results/{{scenario}}/result_summaries/{result_summary}.csv', 
            result_summary = result_summaries),
    log:
        log = 'results/{scenario}/logs/summarise_results.log'
    shell: 
        'python workflow/scripts/osemosys_global/summarise_results.py 2> {log}'