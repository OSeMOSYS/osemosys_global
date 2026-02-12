
rule create_base_generator_data:
    message: "Creating generator table from PLEXOS"
    params:
        start_year = config["startYear"],
    input:
        plexos = "resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx",
        nodes = "resources/data/nodes.yaml", 
        techs = "resources/data/techs.yaml"
    output:
        generators = "results/data/powerplants/generator_template.csv" 
    script:
        "../scripts/osemosys_global/techs/create_generator_template.py"

# rule create_thermal_activity_ratios:
#     message: "Creating IAR/OAR for thermal powerplants"
#     params:
#     input:
#     output:
#     script:


# rule create_residual_capacity:
#     message: "Generating Residual Capacity..."
#     input:
#         plexos = 'resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx',
#         weo_costs = 'resources/data/weo_2020_powerplant_costs.csv',
#         weo_regions = 'resources/data/weo_region_mapping.csv',
#         default_op_life = 'resources/data/operational_life.csv',
#         naming_convention_tech = 'resources/data/naming_convention_tech.csv',
#         line_data = 'resources/data/Costs Line expansion.xlsx',
#         default_av_factors = 'resources/data/availability_factors.csv',
#         custom_res_cap = powerplant_cap_custom_csv()
#     params:
#         trade = config['crossborderTrade'],
#         start_year = config['startYear'],
#         end_year = config['endYear'],
#         invest_techs = config['no_invest_technologies']
#     output:
#         csv_files = expand('results/data/{output_file}', output_file = power_plant_files)
#     log:
#         log = 'results/logs/powerplant.log'
#     script:
#         "../scripts/osemosys_global/powerplant/main.py"