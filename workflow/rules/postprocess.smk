import os

# OUTPUT FILES 

RESULT_FIGURES = [
    'TotalCapacityAnnual', 
    'GenerationAnnual',
]

RESULT_SUMMARIES = [
    "TradeFlows",
    "AnnualEmissionIntensity",
    "TransmissionCapacity",
    "PowerCapacity",
    "GenerationShares",
    "NodeCost"
]

rule otoole_results:
    message:
        'Generating result csv files...'
    input:
        solution_file = "results/{scenario}/{scenario}.sol",
        otoole_config = 'results/{scenario}/otoole.yaml',
    params:
        csv_dir = 'results/{scenario}/data',
    output:
        expand('results/{{scenario}}/results/{result_file}.csv', result_file = OTOOLE_RESULTS),
    log:
        log = 'results/{scenario}/logs/otoole_results.log'
    shell: 
        """
        otoole results {config[solver]} csv \
        {input.solution_file} results/{wildcards.scenario}/results \
        csv {params.csv_dir} {input.otoole_config} 2> {log} 
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
        expand('results/{{scenario}}/figures/{result_figure}.html', result_figure = RESULT_FIGURES)
    log:
        log = 'results/{scenario}/logs/visualisation.log'
    shell: 
        'python workflow/scripts/osemosys_global/visualise.py {params.input_data} {params.result_data} {params.scenario_figs_dir} {params.cost_line_expansion_xlsx} {params.countries} {params.results_by_country} {params.years} 2> {log}'

# rule summarise_results:
#     message:
#         'Generating summary of results...'
#     input:
#         csv_files = expand('results/{{scenario}}/results/{result_file}.csv', result_file = OTOOLE_RESULTS),
#     params:
#         start_year = config['startYear'],
#         end_year = config['endYear'],
#         dayparts = config['dayparts'],
#         seasons = config['seasons'],
#     output:
#         expand('results/{{scenario}}/result_summaries/{result_summary}.csv', 
#             result_summary = result_summaries),
#     log:
#         log = 'results/{scenario}/logs/summarise_results.log'
#     shell: 
#         'python workflow/scripts/osemosys_global/summarise_results.py 2> {log}'

rule calculate_trade_flows:
    message: 
        "Calculating Hourly Trade Flows"
    params:
        seasons = config["seasons"],
        dayparts = config["dayparts"],
        timeshift = config["timeshift"],
    input:
        activity_by_mode = "results/{scenario}/results/TotalAnnualTechnologyActivityByMode.csv",
    output:
        trade_flows = "results/{scenario}/result_summaries/TradeFlows.csv",
    log:
        log = "results/{scenario}/logs/trade_flows.log"
    script: 
        "../scripts/osemosys_global/summary/trade_flows.py"

rule calculate_carbon_intensity:
    message:
        "Calculating Carbon Intensity..."
    input:
        production_by_technology = "results/{scenario}/results/ProductionByTechnologyAnnual.csv",
        annual_emissions = "results/{scenario}/results/AnnualEmissions.csv",
    output:
        emission_intensity = "results/{scenario}/result_summaries/AnnualEmissionIntensity.csv",
    log:
        log = 'results/{scenario}/logs/carbon_intensity.log'
    script: 
        "../scripts/osemosys_global/summary/carbon_intensity.py"

rule calculate_cost_by_node:
    message:
        "Calculating Cost by Node..."
    input:
        discounted_cost_by_technology = "results/{scenario}/results/DiscountedCostByTechnology.csv",
        demand = "results/{scenario}/results/Demand.csv",
    output:
        node_cost = "results/{scenario}/result_summaries/NodeCost.csv",
    log:
        log = 'results/{scenario}/logs/node_cost.log'
    script: 
        "../scripts/osemosys_global/summary/node_cost.py"

rule calculate_generation_shares:
    message:
        "Calculating Generaion Fuel Shares..."
    input:
        production_by_technology = "results/{scenario}/results/ProductionByTechnology.csv",
    output:
        generation_shares = "results/{scenario}/result_summaries/GenerationShares.csv",
    log:
        log = 'results/{scenario}/logs/generation_shares.log'
    script: 
        "../scripts/osemosys_global/summary/gen_shares.py"

rule calculate_capacity_by_node:
    message:
        "Calculating Capacity by Node..."
    input:
        total_capacity = "results/{scenario}/results/TotalCapacityAnnual.csv",
    output:
        power_capacity = "results/{scenario}/result_summaries/PowerCapacity.csv",
        transmission_capacity = "results/{scenario}/result_summaries/TransmissionCapacity.csv",
    log:
        log = 'results/{scenario}/logs/generation_shares.log'
    script: 
        "../scripts/osemosys_global/summary/capacity.py"