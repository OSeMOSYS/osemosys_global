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
    'AvailabilityFactor',
    'TotalAnnualMaxCapacity',
    'AccumulatedAnnualDemand',
    'TotalAnnualMinCapacity'
    ]

transmission_files = [
    'transmission/CapitalCost',
    'transmission/VariableCost',
    'transmission/FixedCost',
    'transmission/CapacityToActivityUnit',
    'transmission/OperationalLife',
    'transmission/TotalAnnualMaxCapacityInvestment',
    'transmission/TotalAnnualMinCapacityInvestment',
    'TotalTechnologyModelPeriodActivityUpperLimit',
    'transmission/InputActivityRatio',
    'transmission/OutputActivityRatio',
    'transmission/ResidualCapacity',
    'transmission/TECHNOLOGY',
    'FUEL'
    ]
    
storage_files = [
    'CapitalCost',
    'CapitalCostStorage',
    'FixedCost',
    'VariableCost',
    'CapacityToActivityUnit',
    'OperationalLife',    
    'OperationalLifeStorage',
    'TotalAnnualMaxCapacityInvestment',
    'TotalAnnualMinCapacityInvestment',
    'InputActivityRatio',
    'OutputActivityRatio',
    'ResidualCapacity',
    'ResidualStorageCapacity',
    'TECHNOLOGY',
    'STORAGE',
    'StorageLevelStart',
    'TechnologyToStorage',
    'TechnologyFromStorage'
    ]

timeslice_files = [
    'CapacityFactor',
    'TIMESLICE',
    'SpecifiedDemandProfile',
    'YearSplit',
    'Conversionls',
    'Conversionld',
    'Conversionlh',
    'SEASON',
    'DAYTYPE',
    'DAILYTIMEBRACKET',
    'DaySplit',
    ]

reserves_files = [
    'ReserveMargin',
    'ReserveMarginTagTechnology',
    'ReserveMarginTagFuel'
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

user_capacity_files = [
    'TotalAnnualMinCapacityInvestment',
    'TotalAnnualMaxCapacityInvestment'
]

fuel_limit_files = [
    'TotalTechnologyAnnualActivityUpperLimit',
]

GENERATED_CSVS = (
    power_plant_files + transmission_files + storage_files + timeslice_files \
    + reserves_files + demand_files + emission_files + fuel_limit_files

)
GENERATED_CSVS = [Path(x).stem for x in GENERATED_CSVS]
EMPTY_CSVS = [x for x in OTOOLE_PARAMS if x not in GENERATED_CSVS]

# rules

rule make_data_dir:
    output: directory('results/data')
    shell: 'mkdir -p {output}'

rule demand_projections:
    message:
        "Generating demand data..."
    input:
        plexos = "resources/data/default/PLEXOS_World_2015_Gold_V1.1.xlsx",
        plexos_demand = "resources/data/default/All_Demand_UTC_2015.csv",
        iamc_gdp ="resources/data/default/iamc_db_GDPppp_Countries.xlsx",
        iamc_pop = "resources/data/default/iamc_db_POP_Countries.xlsx",
        iamc_urb = "resources/data/default/iamc_db_URB_Countries.xlsx",
        iamc_missing = "resources/data/default/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx",
        td_losses = "resources/data/default/T&D Losses.xlsx",
        ember = "resources/data/default/ember_yearly_electricity_data.csv",
        custom_nodes = "resources/data/custom/specified_annual_demand.csv"
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
        
rule powerplant:
    message:
        "Generating powerplant data..."
    input:
        rules.demand_projections.output.csv_files,
        plexos = 'resources/data/default/PLEXOS_World_2015_Gold_V1.1.xlsx',
        res_limit = 'resources/data/default/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx',
        fuel_limit = 'resources/data/custom/fuel_limits.csv',
        build_rates = 'resources/data/custom/powerplant_build_rates.csv',        
        weo_costs = 'resources/data/default/weo_2020_powerplant_costs.csv',
        weo_regions = 'resources/data/default/weo_region_mapping.csv',
        default_op_life = 'resources/data/custom/operational_life.csv',
        naming_convention_tech = 'resources/data/default/naming_convention_tech.csv',
        default_af_factors = 'resources/data/custom/availability_factors.csv',
        custom_res_cap = 'resources/data/custom/residual_capacity.csv',
        custom_res_potentials = 'resources/data/custom/RE_potentials.csv'
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        geographic_scope = config['geographic_scope'],
        custom_nodes = config['nodes_to_add'],
        remove_nodes = config['nodes_to_remove'],
        user_defined_capacity = config['user_defined_capacity'],
        no_investment_techs = config['no_invest_technologies'],
        availability_factors = config['max_availability_factors'],
        res_targets = config['re_targets'],
        fossil_capacity_targets = config['fossil_capacity_targets'],
        calibration = config['min_generation_factors'],
        output_data_dir = 'results/data',
        input_data_dir = 'resources/data/default',
        powerplant_data_dir = 'results/data/powerplant',
    output:
        csv_files = expand('results/data/{output_file}.csv', output_file = power_plant_files)
    log:
        log = 'results/logs/powerplant.log'
    script:
        "../scripts/osemosys_global/powerplant/main.py"

rule powerplant_var_costs:
    message:
        "Generating powerplant variable costs..."
    input:
        cmo_forecasts = 'resources/data/default/CMO-October-2024-Forecasts.xlsx',
        fuel_prices = 'resources/data/custom/fuel_prices.csv',
        regions = "results/data/REGION.csv",
        years = "results/data/YEAR.csv",
        technologies = "results/data/powerplant/TECHNOLOGY.csv",
    output:
        var_costs = 'results/data/powerplant/VariableCost.csv'
    log:
        log = 'results/logs/powerplant_var_cost.log'
    script:
        "../scripts/osemosys_global/powerplant/variable_costs.py"

rule fuel_limits:
    message:
        "Generating mining fuel limits..."
    input:
        region_csv = "results/data/REGION.csv",
        technology_csv = "results/data/powerplant/TECHNOLOGY.csv",
        year_csv = "results/data/YEAR.csv",
        fuel_limit_csv = "resources/data/custom/fuel_limits.csv",
    output:
        activity_upper_limit_csv = 'results/data/TotalTechnologyAnnualActivityUpperLimit.csv'
    log:
        log = 'results/logs/powerplant_fuel_limits.log'
    script:
        "../scripts/osemosys_global/powerplant/fuel_limits.py"

rule transmission:
    message:
        "Generating transmission data..."
    input:
        rules.powerplant.output.csv_files,
        'results/data/powerplant/VariableCost.csv',
        default_op_life = 'resources/data/custom/operational_life.csv',
        gtd_existing = 'resources/data/default/GTD_existing.csv',
        gtd_planned = 'resources/data/default/GTD_planned.csv',
        gtd_mapping = 'resources/data/default/GTD_region_mapping.csv',
        centerpoints = 'resources/data/default/centerpoints.csv',
        transmission_build_rates = 'resources/data/custom/transmission_build_rates.csv',
    params:
        trade = config['crossborderTrade'],
        transmission_existing = config['transmission_existing'],
        transmission_planned = config['transmission_planned'],
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        custom_nodes = config['nodes_to_add'],
        user_defined_capacity_transmission = config['user_defined_capacity_transmission'],
        no_investment_techs = config['no_invest_technologies'],
        transmission_parameters = config['transmission_parameters'],
        output_data_dir = 'results/data',
        powerplant_data_dir = 'results/data/powerplant',
        transmission_data_dir = 'results/data/transmission',
    output:
        csv_files = expand('results/data/{output_file}.csv', output_file = transmission_files)
    log:
        log = 'results/logs/transmission.log'
    script:
        "../scripts/osemosys_global/transmission/main.py"
        
rule storage:
    message:
        "Generating storage data..."
    input:
        rules.transmission.output.csv_files,
        default_op_life = 'resources/data/custom/operational_life.csv',
        storage_build_rates = 'resources/data/custom/storage_build_rates.csv',
        gesdb_project_data = 'resources/data/default/GESDB_Project_Data.json',
        gesdb_regional_mapping = 'resources/data/default/GESDB_region_mapping.csv',
    params:
        storage_existing = config['storage_existing'],
        storage_planned = config['storage_planned'],
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        custom_nodes = config['nodes_to_add'],
        user_defined_capacity_storage = config['user_defined_capacity_storage'],
        no_investment_techs = config['no_invest_technologies'],
        storage_parameters = config['storage_parameters'],
        output_data_dir = 'results/data',
        transmission_data_dir = 'results/data/transmission',
    output:
        csv_files = expand('results/data/{output_file}.csv', output_file = storage_files)
    log:
        log = 'results/logs/storage.log'
    script:
        "../scripts/osemosys_global/storage/main.py"        

rule timeslice:
    message:
        'Generating timeslice data...'
    input:
        plexos_demand = 'resources/data/default/All_Demand_UTC_2015.csv',
        plexos_csp_2015 = 'resources/data/default/CSP 2015.csv',
        plexos_spv_2015 = 'resources/data/default/SolarPV 2015.csv',
        plexos_hyd_2015 = 'resources/data/default/Hydro_Monthly_Profiles (15 year average).csv',
        plexos_won_2015 = 'resources/data/default/Won 2015.csv',
        plexos_wof_2015 = 'resources/data/default/Woff 2015.csv',
        custom_specified_demand_profiles = 'resources/data/custom/specified_demand_profile.csv',
        custom_csp_profiles = 'resources/data/custom/RE_profiles_CSP.csv',
        custom_hyd_profiles = 'resources/data/custom/RE_profiles_HYD.csv',
        custom_spv_profiles = 'resources/data/custom/RE_profiles_SPV.csv',
        custom_wof_profiles = 'resources/data/custom/RE_profiles_WOF.csv',
        custom_won_profiles = 'resources/data/custom/RE_profiles_WON.csv',                                
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        output_data_dir = 'results/data',
        input_data_dir = 'resources/data/default',
        input_dir = 'resources',
        output_dir = 'results',
        custom_nodes_dir = 'resources/data/custom',
        geographic_scope = config['geographic_scope'],
        seasons = config['seasons'],
        dayparts = config['dayparts'],
        daytype = config['daytype'],        
        timeshift = config['timeshift'],
    output:
        csv_files = expand('results/data/{output_file}.csv', output_file=timeslice_files),
    log:
        log = 'results/logs/timeslice.log'    
    script:
        "../scripts/osemosys_global/TS_data.py"
        
rule reserves:
    message:
        'Generating reserves data...'
    input:
        rules.storage.output.csv_files,
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        output_data_dir = 'results/data',
        reserve_margin = config['reserve_margin'],
        reserve_margin_technologies = config['reserve_margin_technologies']
    output:
        csv_files = expand('results/data/{output_file}.csv', output_file=reserves_files),
    log:
        log = 'results/logs/reserves.log'    
    script:
        "../scripts/osemosys_global/reserves/main.py"       

rule demand_projection_figures:
    message:
        "Generating demand figures..."
    input:
        plexos = "resources/data/default/PLEXOS_World_2015_Gold_V1.1.xlsx",
        iamc_gdp ="resources/data/default/iamc_db_GDPppp_Countries.xlsx",
        iamc_pop = "resources/data/default/iamc_db_POP_Countries.xlsx",
        iamc_urb = "resources/data/default/iamc_db_URB_Countries.xlsx",
        iamc_missing = "resources/data/default/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx",
        ember = "resources/data/default/ember_yearly_electricity_data.csv"
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
        ember = 'resources/data/default/ember_yearly_electricity_data.csv',    
        emissions_factors = 'resources/data/default/emission_factors.csv',
        iar = 'results/data/InputActivityRatio.csv',
        oar = 'results/data/OutputActivityRatio.csv',
    params:
        start_year = config['startYear'],
        end_year = config['endYear'],
        region_name = 'GLOBAL',
        output_data_dir = 'results/data',
        emission_penalty = config['emission_penalty'],
        emission_limit = config['emission_limit'],
    output: 
        csv_files = expand('results/data/{output_file}.csv', output_file = emission_files),
    log:
        log = 'results/logs/emissions.log'
    script:
        "../scripts/osemosys_global/emissions/main.py"   

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