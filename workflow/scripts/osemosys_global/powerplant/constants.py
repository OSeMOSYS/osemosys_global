"""Constants for the powerplant module"""

"""Change IAR for CSP value taken from PLEXOS to 1.0. 
Relevant for the 'average_efficiency' function."""
avg_csp_eff = 1

"""Change IAR for URN value taken from PLEXOS to 2.2 (45%). 
Relevant for the 'average_efficiency' function."""
avg_urn_eff = 0.45

"""Technologies that will have 00 and 01 suffixes to represent PLEXOS 
   historical values and future values. Relevant for the 'residual_capacity' 
   and activity functions."""
duplicate_techs = ['CCG', 'OCG']
    
"""Add extra nodes which exist in 2050 but are not in the 2015 PLEXOS-World 
data to enable their addition to the workflow. Relevant for the 
'generator_table' and activity functions"""
nodes_extra_list = ['AF-SOM',
                    'AF-TCD',
                    'AS-TLS',
                    'EU-MLT',
                    'NA-BLZ',
                    'NA-HTI',
                    'SA-BRA-J1',
                    'SA-BRA-J2',
                    'SA-BRA-J3',
                    'SA-SUR']

"""Sets the mode of operations for technologies. 
Relevant for the activity_master_start function."""
mode_list = [1,2]

"""List of technologies (thermal) for which output activity ratios need 
to be developed. Relevant for the 'activity_output_pwr' function."""
thermal_fuel_list_oar = ['COA',
                     'COG',
                     'OCG',
                     'CCG',
                     'PET',
                     'URN',
                     'OIL',
                     'OTH',
                     'CCS'
                    ]

"""List of technologies (thermal) for which input activity ratios need 
to be developed. Relevant for the 'activity_input_pwr' function."""
thermal_fuel_list_iar = ['COA',
                         'COG',
                         'PET',
                         'URN',
                         'OIL',
                         'OTH'
                        ]

"""List of upstream fuels for which output activity ratios need 
to be developed. Relevant for the 'activity_upstream' function."""
thermal_fuel_list_mining = ['COA',
                             'COG',
                             'GAS',
                             'PET',
                             'URN',
                             'OIL',
                             'OTH'
                            ]

"""List of technologies (renewable) for which activity ratios need 
to be developed. Relevant for the 'activity_input_pwr' and 
'activity_upstream' functions."""
renewables_list = ['BIO',
                   'GEO',
                   'HYD',
                   'SPV', 
                   'CSP',
                   'WAS',
                   'WAV',
                   'WON', 
                   'WOF']

"""Technology input for the mapping of WEO cost data to OG technologies.
Relevant for the 'costs_pwr' function."""
costs_dict = {'Biomass - waste incineration - CHP':'WAS',
              'Biomass Power plant':'BIO', 
              'CCGT':'CCG', 
              'CCGT - CHP':'COG', 
              'Concentrating solar power':'CSP',
              'Gas turbine':'OCG',
              'Geothermal':'GEO', 
              'Hydropower - large-scale':'HYD',
              'Marine':'WAV',
              'Nuclear':'URN', 
              'Solar photovoltaics - Large scale':'SPV', 
              'Steam Coal - SUBCRITICAL':'COA',
              'Steam Coal - SUPERCRITICAL':'COA', 
              'Steam Coal - ULTRASUPERCRITICAL':'COA',
              'Wind onshore':'WON',
              'Wind offshore':'WOF',
              'Petroleum':'PET',
              'Oil':'OIL',
              'Other':'OTH',
              'IGCC + CCS':'CCS',
              'Coal* + CCS':'CCS',}

"""setings for custom iar values deviating from derived values from the 
   PLEXOS-World dataset. Relevant for the 'newIar' function."""
new_iar_ccg = 0.5
new_iar_ocg = 0.35
new_iar_coa = 0.33
new_iar_default = 1

"""Set iar and oar values for custom transmission entries. I.e. oar of 0.9
assumes 10% losses. Relevant for the user_defined_capacity function."""
df_iar_custom_val = 1
df_oar_custom_val = 0.9

"""Set column name dictionaries for different Global Energy Monitor (gem) input datasets"""
gem_coal_col = {'Country' : 'Country', 'Capacity (MW)' : 'VALUE', 
                'Status' : 'Status', 'Year' : 'Year_built',
                'RETIRED' : 'Year_retired', 'Planned Retire' : 'Year_retired_planned', 
                'Latitude' : 'Latitude', 'Longitude' : 'Longitude'}

gem_gas_col = {'Country' : 'Country', 'Capacity elec. (MW)' : 'VALUE', 
                'Status' : 'Status', 'Start year' : 'Year_built',
                'Retired year' : 'Year_retired', 'Planned retire' : 'Year_retired_planned', 
                'Latitude' : 'Latitude', 'Longitude' : 'Longitude', 
                'Technology' : 'Technology'}

"""Set technology dictionary to match with OSeMOSYS global technologies"""
gem_gas_dict = {'CC' : 'CCG', 
                'GT' : 'OCG',
                'ICCC' : 'CCG',
                'ISCC' : 'CCG',
                'ST' : 'OCG',
                'AFC' : 'CCG'}

if "snakemake" in globals():
    file_plexos = snakemake.input.plexos
    file_default_op_life = snakemake.input.default_op_life
    file_naming_convention_tech = snakemake.input.naming_convention_tech
    file_weo_costs = snakemake.input.weo_costs
    file_weo_regions = snakemake.input.weo_regions
    file_line_data = snakemake.input.line_data
    file_default_av_factors = snakemake.input.default_av_factors
    file_custom_res_cap = snakemake.input.custom_res_cap
    start_year = snakemake.params.start_year
    end_year = snakemake.params.end_year
    region_name = snakemake.params.region_name
    custom_nodes = snakemake.params.custom_nodes
    custom_nodes_data = snakemake.input.custom_nodes
    tech_capacity = snakemake.params.user_defined_capacity
    cross_border_trade = snakemake.params.crossborderTrade
    no_investment_techs = snakemake.params.no_invest_technologies
    output_data_dir = snakemake.params.output_data_dir
    input_data_dir = snakemake.params.output_data_dir

else:
    file_plexos = 'resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx'
    file_default_op_life = 'resources/data/operational_life.csv'
    file_naming_convention_tech = 'resources/data/naming_convention_tech.csv'
    file_weo_costs = 'resources/data/weo_2020_powerplant_costs.csv'
    file_weo_regions = 'resources/data/weo_region_mapping.csv'
    file_line_data = 'resources/data/Costs Line expansion.xlsx'
    file_default_av_factors = 'resources/data/availability_factors.csv'    
    file_custom_res_cap = 'resources/data/custom_nodes/residual_capacity.csv'     
    start_year = 2020
    end_year = 2050
    region_name = 'GLOBAL'
    custom_nodes = ["INDWE", "INDEA", "INDNE", "INDNO", "INDSO"]
    custom_nodes_data = 'resources/data/custom_nodes/specified_annual_demand.csv'
    tech_capacity = {'TRNINDEAINDNE': [0, 2030, "open", 2030, 10, 861]}
    cross_border_trade = True
    no_investment_techs = ["CSP", "WAV", "URN", "OTH", "WAS", 
                           "COG", "GEO", "BIO", "PET"]
    output_data_dir = 'results/data'
    input_data_dir = 'resources/data'
    
years = list(range(start_year, end_year + 1))

SET_DTYPES = {
    "DAILYTIMEBRACKET": int,
    "EMISSION":str,
    "FUEL":str,
    "MODE_OF_OPERATION":int,
    "REGION":str,
    "SEASON":str,
    "STORAGE":str,
    "TECHNOLOGY":str,
    "TIMESLICE":str,
    "YEAR":int,
}