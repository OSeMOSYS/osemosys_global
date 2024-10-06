import os

# OUTPUT FILES 

result_figures = [
    'TotalCapacityAnnual', 
    'GenerationAnnual',
]

result_summaries = [
    'Metrics'
]

# rules

rule otoole_results:
    message:
        'Generating result csv files...'
    input:
        solution_file = "results/{scenario}/{scenario}.sol",
        datafile = 'results/{scenario}/{scenario}.txt',
        otoole_config = 'results/{scenario}/otoole.yaml',
    output:
        expand('results/{{scenario}}/results/{result_file}.csv', result_file = OTOOLE_RESULTS),
    # conda:
    #     '../envs/otoole.yaml'
    log:
        log = 'results/{scenario}/logs/otoole_results.log'
    shell: 
        """
        otoole results {config[solver]} csv \
        {input.solution_file} results/{wildcards.scenario}/results \
        datafile {input.datafile} \
        {input.otoole_config}
        2> {log} 
        """

rule visualisation:
    message:
        'Generating result figures...'
    input:
        csv_files = expand('results/{{scenario}}/results/{result_file}.csv', result_file = OTOOLE_RESULTS),
    params:
        input_data = "results/{scenario}/data/",
        result_data = "results/{scenario}/results/",
        scenario_figs_dir = "results/{scenario}/figures/",
        cost_line_expansion_xlsx = "'resources/data/Costs Line expansion.xlsx'",
        countries = config['geographic_scope'],
        results_by_country = config['results_by_country'],
        years = [config['endYear']],
    output:
        expand('results/{{scenario}}/figures/{result_figure}.html', result_figure = result_figures)
    log:
        log = 'results/{scenario}/logs/visualisation.log'
    shell: 
        'python workflow/scripts/osemosys_global/visualise.py {params.input_data} {params.result_data} {params.scenario_figs_dir} {params.cost_line_expansion_xlsx} {params.countries} {params.results_by_country} {params.years} 2> {log}'

rule summarise_results:
    message:
        'Generating summary of results...'
    input:
        csv_files = expand('results/{{scenario}}/results/{result_file}.csv', result_file = OTOOLE_RESULTS),
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