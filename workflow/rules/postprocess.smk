import os

# OUTPUT FILES 

RESULT_FIGURES = [
    'TotalCapacityAnnual', 
    'GenerationAnnual',
]

RESULT_SUMMARIES = [
    "TradeFlowsNode",
    "TradeFlowsCountry",
    "AnnualNetTradeFlowsNode",
    "AnnualNetTradeFlowsCountry",
    "AnnualImportTradeFlowsNode",
    "AnnualImportTradeFlowsCountry",
    "AnnualExportTradeFlowsNode",
    "AnnualExportTradeFlowsCountry",
    "AnnualTotalTradeFlowsNode",
    "AnnualTotalTradeFlowsCountry",
    "AnnualEmissionIntensity",
    "PowerCapacityNode",
    "TransmissionCapacityNode",
    "PowerCapacityCountry",
    "TransmissionCapacityCountry",
    "GenerationSharesNode",
    "GenerationSharesCountry",
    "GenerationSharesGlobal",
    "PowerCostNode",
    "TotalCostNode",
    "PowerCostCountry",
    "TotalCostCountry",
    "PowerCostGlobal",
    "TotalCostGlobal",
    "Metrics"
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
        centerpoints = 'resources/data/centerpoints.csv',
        custom_nodes_centerpoints = 'resources/data/custom_nodes/centerpoints.csv',
        color_codes = 'resources/data/color_codes.csv',
    params:
        result_input_data = "results/{scenario}/data/",
        result_data = "results/{scenario}/results/",
        scenario_figs_dir = "results/{scenario}/figures/",
        geographic_scope = config['geographic_scope'],
        results_by_country = config['results_by_country'],
        start_year = config['startYear'],
        end_year = [config['endYear']],
        custom_nodes = config['nodes_to_add'],
        seasons = config['seasons'],
        dayparts = config['dayparts'],
        timeshift = config['timeshift'],
    output:
        expand('results/{{scenario}}/figures/{result_figure}.html', result_figure = RESULT_FIGURES)
    log:
        log = 'results/{scenario}/logs/visualisation.log'
    script: 
        "../scripts/osemosys_global/visualisation/visualise.py" 

rule calculate_trade_flows:
    message: 
        "Calculating Trade Flows..."
    params:
        seasons = config["seasons"],
        dayparts = config["dayparts"],
        timeshift = config["timeshift"],
    input:
        activity_by_mode = "results/{scenario}/results/TotalAnnualTechnologyActivityByMode.csv",
    output:
        node_trade_flows = "results/{scenario}/result_summaries/TradeFlowsNode.csv",
        country_trade_flows = "results/{scenario}/result_summaries/TradeFlowsCountry.csv",
        annual_net_node_trade_flows = "results/{scenario}/result_summaries/AnnualNetTradeFlowsNode.csv",
        annual_net_country_trade_flows = "results/{scenario}/result_summaries/AnnualNetTradeFlowsCountry.csv",
        annual_import_node_trade_flows = "results/{scenario}/result_summaries/AnnualImportTradeFlowsNode.csv",
        annual_import_country_trade_flows = "results/{scenario}/result_summaries/AnnualImportTradeFlowsCountry.csv",        
        annual_export_node_trade_flows = "results/{scenario}/result_summaries/AnnualExportTradeFlowsNode.csv",
        annual_export_country_trade_flows = "results/{scenario}/result_summaries/AnnualExportTradeFlowsCountry.csv",        
        annual_total_node_trade_flows = "results/{scenario}/result_summaries/AnnualTotalTradeFlowsNode.csv",
        annual_total_country_trade_flows = "results/{scenario}/result_summaries/AnnualTotalTradeFlowsCountry.csv",        
    log:
        log = "results/{scenario}/logs/trade_flows.log"
    script: 
        "../scripts/osemosys_global/summary/trade_flows.py"

rule calculate_carbon_intensity:
    message:
        "Calculating Carbon Intensity..."
    params:
        storage = config['storage_parameters'],
    input:
        production_by_technology = "results/{scenario}/results/ProductionByTechnology.csv",
        annual_emissions = "results/{scenario}/results/AnnualEmissions.csv",
    output:
        emission_intensity = "results/{scenario}/result_summaries/AnnualEmissionIntensity.csv",
    log:
        log = 'results/{scenario}/logs/carbon_intensity.log'
    script: 
        "../scripts/osemosys_global/summary/carbon_intensity.py"

rule calculate_costs:
    message:
        "Calculating Costs..."
    input:
        discounted_cost_by_technology = "results/{scenario}/results/DiscountedCostByTechnology.csv",
        demand = "results/{scenario}/results/Demand.csv",
    output:
        node_pwr_cost = "results/{scenario}/result_summaries/PowerCostNode.csv",
        country_pwr_cost = "results/{scenario}/result_summaries/PowerCostCountry.csv",
        global_pwr_cost = "results/{scenario}/result_summaries/PowerCostGlobal.csv",
        node_cost = "results/{scenario}/result_summaries/TotalCostNode.csv",
        country_cost = "results/{scenario}/result_summaries/TotalCostCountry.csv",
        global_cost = "results/{scenario}/result_summaries/TotalCostGlobal.csv",
    log:
        log = 'results/{scenario}/logs/node_cost.log'
    script: 
        "../scripts/osemosys_global/summary/costs.py"

rule calculate_generation_shares:
    message:
        "Calculating Generaion Fuel Shares..."
    params:
        storage = config['storage_parameters'],
    input:
        production_by_technology = "results/{scenario}/results/ProductionByTechnology.csv",
    output:
        generation_shares_node = "results/{scenario}/result_summaries/GenerationSharesNode.csv",
        generation_shares_country = "results/{scenario}/result_summaries/GenerationSharesCountry.csv",
        generation_shares_global = "results/{scenario}/result_summaries/GenerationSharesGlobal.csv",
    log:
        log = 'results/{scenario}/logs/generation_shares.log'
    script: 
        "../scripts/osemosys_global/summary/gen_shares.py"

rule calculate_capacity:
    message:
        "Calculating Capacities..."
    input:
        total_capacity = "results/{scenario}/results/TotalCapacityAnnual.csv",
    output:
        power_capacity_node = "results/{scenario}/result_summaries/PowerCapacityNode.csv",
        transmission_capacity_node = "results/{scenario}/result_summaries/TransmissionCapacityNode.csv",
        power_capacity_country = "results/{scenario}/result_summaries/PowerCapacityCountry.csv",
        transmission_capacity_country = "results/{scenario}/result_summaries/TransmissionCapacityCountry.csv",
    log:
        log = 'results/{scenario}/logs/generation_shares.log'
    script: 
        "../scripts/osemosys_global/summary/capacity.py"

rule calcualte_headline_metrics:
    message:
        "Calculating Headline Metrics..."
    params:
        storage = config['storage_parameters'],
    input:
        annual_emissions = "results/{scenario}/results/AnnualEmissions.csv",
        production_by_technology = "results/{scenario}/results/ProductionByTechnologyAnnual.csv",
        total_discounted_cost = "results/{scenario}/results/TotalDiscountedCost.csv",
        demand = "results/{scenario}/results/Demand.csv"
    output:
        metrics = "results/{scenario}/result_summaries/Metrics.csv",
    log:
        log = 'results/{scenario}/logs/generation_shares.log'
    script: 
        "../scripts/osemosys_global/summary/headline.py"