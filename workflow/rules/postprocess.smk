"""Rules for processing OSeMOSYS Solutions"""

import os
configfile: 'config/config.yaml'

# input functions

def solver_file_type(wildcards):
    if config['solver'] == 'cplex':
        return 'results/{scenario}/{scenario}_sort.sol'
    else: 
        return 'results/{scenario}/{scenario}.sol'

# constants 

RESULT_FIGURES = [
    'TotalCapacityAnnual', 
    'GenerationAnnual',
]

RESULT_SUMMARIES = [
    'Metrics'
]

OUTPUT_FILES = pd.read_csv('resources/otoole_files.csv')['outputs'].dropna().to_list()

# rules

rule otoole_results:
    message:
        'Generating result csv files...'
    input:
        solution_file = solver_file_type,
        pre_process_file = 'results/{scenario}/PreProcessed_{scenario}.txt'
    params: 
        datapackage = 'results/{scenario}/datapackage.json',
        config = 'results/{scenario}/config.yaml',
    output:
        expand('results/{{scenario}}/results/{result_file}.csv', result_file = OUTPUT_FILES),
    conda:
        '../envs/otoole.yaml'
    log:
        log = 'results/{scenario}/logs/otoole_results.log'
    shell: 
        '''
        otoole results {config[solver]} csv \
        {input.solution_file} results/{wildcards.scenario}/results \
        --input_datafile {input.pre_process_file} \
        --input_datapackage {params.datapackage} \
        {params.config} 2> {log} 
        '''

rule visualisation:
    message:
        'Generating result figures...'
    input:
        csv_files = expand('results/{{scenario}}/results/{result_file}.csv', result_file = OUTPUT_FILES),
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        dayparts = config['dayparts'],
        seasons = config['seasons'],
        by_country = config['results_by_country'],
        geographic_scope = config['geographic_scope'],
    output:
        expand('results/{{scenario}}/figures/{result_figure}.html', result_figure = RESULT_FIGURES)
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/{scenario}/logs/visualisation.log'
    shell: 
        'python workflow/scripts/osemosys_global/visualisation.py 2> {log}'

rule summarise_results:
    message:
        'Generating summary of results...'
    input:
        csv_files = expand('results/{{scenario}}/results/{result_file}.csv', result_file = OUTPUT_FILES),
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        dayparts = config['dayparts'],
        seasons = config['seasons'],
    output:
        expand('results/{{scenario}}/result_summaries/{result_summary}.csv', 
            result_summary = RESULT_SUMMARIES),
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/{scenario}/logs/summarise_results.log'
    shell: 
        'python workflow/scripts/osemosys_global/summarise_results.py 2> {log}'