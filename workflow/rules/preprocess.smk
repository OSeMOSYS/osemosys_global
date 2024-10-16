import os
import yaml

# model and scenario output files 

demand_figures = [
    'South America',
    'Oceania',
    'North America',
    'Europe',
    'Asia',
    'Africa'
]

# output script files

power_plant_files = [
    'powerplant/CapitalCost',
    'powerplant/FixedCost',
    'powerplant/CapacityToActivityUnit',
    'powerplant/OperationalLife',
    'powerplant/TotalAnnualMaxCapacityInvestment',
    'powerplant/TotalAnnualMinCapacityInvestment',
    'powerplant/FUEL',
    'powerplant/InputActivityRatio',
    'powerplant/OutputActivityRatio',
    'MODE_OF_OPERATION',
    'REGION',
    'powerplant/ResidualCapacity',
    'powerplant/TECHNOLOGY',
    'YEAR',
    'AvailabilityFactor'
    ]

transmission_files = [
    'CapitalCost',
    'FixedCost',
    'CapacityToActivityUnit',
    'OperationalLife',
    'TotalAnnualMaxCapacityInvestment',
    'TotalAnnualMinCapacityInvestment',
    'TotalTechnologyModelPeriodActivityUpperLimit',
    'InputActivityRatio',
    'OutputActivityRatio',
    'ResidualCapacity',
    'TECHNOLOGY',
    'FUEL'
    ]

timeslice_files = [
    'CapacityFactor',
    'TIMESLICE',
    'SpecifiedDemandProfile',
    'YearSplit',
    'STORAGE',
    'TechnologyToStorage',
    'TechnologyFromStorage',
    'Conversionls',
    'Conversionld',
    'Conversionlh',
    'SEASON',
    'DAYTYPE',
    'DAILYTIMEBRACKET',
    'CapitalCostStorage',
    'DaySplit',
    'ReserveMargin',
    'ReserveMarginTagTechnology',
    'ReserveMarginTagFuel'
    ]

variable_cost_files = [
    'VariableCost'
    ]

demand_files = [
    'SpecifiedAnnualDemand'
    ]

emission_files = [
    'EmissionActivityRatio',
    'EmissionsPenalty',
    'EMISSION',
    'AnnualEmissionLimit'
]

max_capacity_files = [
    'TotalAnnualMaxCapacity',
    'TotalTechnologyAnnualActivityUpperLimit',
    'AccumulatedAnnualDemand',
#    'TotalTechnologyModelPeriodActivityUpperLimit'
]

user_capacity_files = [
    'TotalAnnualMinCapacityInvestment',
    'TotalAnnualMaxCapacityInvestment'
]

GENERATED_CSVS = (
    power_plant_files + transmission_files + timeslice_files + variable_cost_files \
    + demand_files + emission_files + max_capacity_files
)
GENERATED_CSVS = [Path(x).stem for x in GENERATED_CSVS]
EMPTY_CSVS = [x for x in OTOOLE_PARAMS if x not in GENERATED_CSVS]

# rules

rule make_data_dir:
    output: directory('results/data')
    shell: 'mkdir -p {output}'

     
def powerplant_cap_custom_csv() -> str:
    if config["nodes_to_add"]:
        return "resources/data/custom_nodes/residual_capacity.csv"
    else:
        return []   
        
rule powerplant:
    message:
        "Generating powerplant data..."
    input:
        plexos = 'resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx',
        weo_costs = 'resources/data/weo_2020_powerplant_costs.csv',
        weo_regions = 'resources/data/weo_region_mapping.csv',
        default_op_life = 'resources/data/operational_life.csv',
        naming_convention_tech = 'resources/data/naming_convention_tech.csv',
        default_av_factors = 'resources/data/availability_factors.csv',
        custom_res_cap = powerplant_cap_custom_csv()
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        custom_nodes = config['nodes_to_add'],
        user_defined_capacity = config['user_defined_capacity'],
        no_investment_techs = config['no_invest_technologies'],
        output_data_dir = 'results/data',
        input_data_dir = 'resources/data',
        powerplant_data_dir = 'results/data/powerplant',

    output:
        csv_files = expand('results/data/{output_file}.csv', output_file = power_plant_files)
    log:
        log = 'results/logs/powerplant.log'
    script:
        "../scripts/osemosys_global/powerplant/main.py"

rule transmission:
    message:
        "Generating transmission data..."
    input:
        rules.powerplant.output.csv_files,
        default_op_life = 'resources/data/operational_life.csv',
        gtd_existing = 'resources/data/GTD_existing.csv',
        gtd_planned = 'resources/data/GTD_planned.csv',
        gtd_mapping = 'resources/data/GTD_region_mapping.csv',
        centerpoints = 'resources/data/centerpoints.csv'
    params:
        trade = config['crossborderTrade'],
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        custom_nodes = config['nodes_to_add'],
        user_defined_capacity_transmission = config['user_defined_capacity_transmission'],
        no_investment_techs = config['no_invest_technologies'],
        transmission_parameters = config['transmission_parameters'],
        output_data_dir = 'results/data',
        input_data_dir = 'resources/data',
        powerplant_data_dir = 'results/data/powerplant',
    output:
        csv_files = expand('results/data/{output_file}.csv', output_file = transmission_files)
    log:
        log = 'results/logs/transmission.log'
    script:
        "../scripts/osemosys_global/transmission/main.py"

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
        csv_files = expand('results/data/{output_file}.csv', output_file=timeslice_files),
    log:
        log = 'results/logs/timeslice.log'    
    shell:
        'python workflow/scripts/osemosys_global/TS_data.py 2> {log}'

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
        csv_files = expand('results/data/{output_file}.csv', output_file=variable_cost_files),
    log:
        log = 'results/logs/variable_costs.log'
    shell:
        'python workflow/scripts/osemosys_global/variablecosts.py 2> {log}'

def demand_custom_csv() -> str:
    if config["nodes_to_add"]:
        return "resources/data/custom_nodes/specified_annual_demand.csv"
    else:
        return []

rule demand_projections:
    message:
        "Generating demand data..."
    input:
        plexos = "resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx",
        plexos_demand = "resources/data/All_Demand_UTC_2015.csv",
        iamc_gdp ="resources/data/iamc_db_GDPppp_Countries.xlsx",
        iamc_pop = "resources/data/iamc_db_POP_Countries.xlsx",
        iamc_urb = "resources/data/iamc_db_URB_Countries.xlsx",
        iamc_missing = "resources/data/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx",
        td_losses = "resources/data/T&D Losses.xlsx",
        ember = "resources/data/ember_yearly_electricity_data.csv",
        custom_nodes = demand_custom_csv()
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        custom_nodes = config["nodes_to_add"]
    output:
        csv_files = 'results/data/SpecifiedAnnualDemand.csv',
    log:
        log = 'results/logs/demand_projections.log'
    script:
        "../scripts/osemosys_global/demand/main.py"

rule demand_projection_figures:
    message:
        "Generating demand figures..."
    input:
        plexos = "resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx",
        iamc_gdp ="resources/data/iamc_db_GDPppp_Countries.xlsx",
        iamc_pop = "resources/data/iamc_db_POP_Countries.xlsx",
        iamc_urb = "resources/data/iamc_db_URB_Countries.xlsx",
        iamc_missing = "resources/data/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx",
        ember = "resources/data/ember_yearly_electricity_data.csv"
    output:
        regression = 'results/figs/regression.png',
        projection = 'results/figs/projection.png'
    log:
        log = 'results/logs/demand_projection_plot.log'
    script:
        "../scripts/osemosys_global/demand/figures.py"

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
        csv_files = expand('results/data/{output_file}.csv', output_file = emission_files),
    log:
        log = 'results/logs/emissions.log'
    shell:
        'python workflow/scripts/osemosys_global/emissions.py 2> {log}'

rule max_capacity:
    message: 
        'Generating capacity limits...'
    input:
        'resources/data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx',
        'results/data/ResidualCapacity.csv',
        'results/data/SpecifiedAnnualDemand.csv'
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
    output:
        csv_files = expand('results/data/{output_file}.csv', output_file = max_capacity_files),
    log:
        log = 'results/logs/max_capacity.log'
    shell:
        'python workflow/scripts/osemosys_global/max_capacity.py 2> {log}'

rule create_missing_csv:
    message:
        "Creating empty parameter data"
    params:
        out_dir = "results/data/"
    input:
        otoole_config = OTOOLE_YAML,
        csvs = expand("results/data/{full}.csv", full=GENERATED_CSVS)
    output:
        csvs = expand("results/data/{empty}.csv", empty=EMPTY_CSVS)
    script:
        "../scripts/osemosys_global/create_missing_csvs.py"
