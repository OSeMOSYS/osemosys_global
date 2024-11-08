"""Constants for the powerplant module"""

"""Change IAR for CSP value taken from PLEXOS to 1.0. 
Relevant for the 'average_efficiency' function."""
AVG_CSP_EFF = 1

"""Change IAR for URN value taken from PLEXOS to 2.2 (45%). 
Relevant for the 'average_efficiency' function."""
AVG_URN_EFF = 0.45

"""Set the year for which data is taken to calculate 
fuel prices/variable costs if to be retrieved from the 
World Bank Commodity Price Forecasts (CMO)."""
CMO_DATA_YEAR = 2023

'''Set default variable cost for techs that do not exist in CMO 
in m$/PJ ($/GJ).'''
BIOMASS_VAR_COSTS = 10
NUCLEAR_VAR_COSTS = 2
WASTE_VAR_COSTS = 10

'''Set the cost factor that scales fuel prices from domestic costs to 
international costs with assumed costs for transportation and shipping.'''
INT_COST_FACTOR = 1.3

"""Technologies that will have 00 and 01 suffixes to represent PLEXOS 
   historical values and future values. Relevant for the 'residual_capacity' 
   and activity functions."""
DUPLICATE_TECHS = ['CCG', 'OCG']
    
"""Add extra nodes which exist in 2050 but are not in the 2015 PLEXOS-World 
data to enable their addition to the workflow. Relevant for the 
'generator_table' and activity functions"""
NODES_EXTRA_LIST = ['AF-SOM',
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
MODE_LIST = [1,2]

"""List of technologies (thermal) for which output activity ratios need 
to be developed. Relevant for the 'activity_output_pwr' function."""
THERMAL_FUEL_LIST_OAR = ['COA',
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
THERMAL_FUEL_LIST_IAR = ['COA',
                         'COG',
                         'PET',
                         'URN',
                         'OIL',
                         'OTH'
                        ]

"""List of upstream fuels for which output activity ratios need 
to be developed. Relevant for the 'activity_upstream' function."""
THERMAL_FUEL_LIST_MINING = ['COA',
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
RENEWABLES_LIST = ['BIO',
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
COSTS_DICT = {'Biomass - waste incineration - CHP':'WAS',
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
NEW_IAR_CCG = 0.5
NEW_IAR_OCG = 0.35
NEW_IAR_COA = 0.33
NEW_IAR_DEFAULT = 1

"""Set column name dictionaries for different Global Energy Monitor (gem) input datasets"""
GEM_COAL_COL = {'Country' : 'Country', 'Capacity (MW)' : 'VALUE', 
                'Status' : 'Status', 'Year' : 'Year_built',
                'RETIRED' : 'Year_retired', 'Planned Retire' : 'Year_retired_planned', 
                'Latitude' : 'Latitude', 'Longitude' : 'Longitude'}

GEM_GAS_COL = {'Country' : 'Country', 'Capacity elec. (MW)' : 'VALUE', 
                'Status' : 'Status', 'Start year' : 'Year_built',
                'Retired year' : 'Year_retired', 'Planned retire' : 'Year_retired_planned', 
                'Latitude' : 'Latitude', 'Longitude' : 'Longitude', 
                'Technology' : 'Technology'}

"""Set technology dictionary to match with OSeMOSYS global technologies"""
GEM_GAS_DICT = {'CC' : 'CCG', 
                'GT' : 'OCG',
                'ICCC' : 'CCG',
                'ISCC' : 'CCG',
                'ST' : 'OCG',
                'AFC' : 'CCG'}

"""Set technology dictionary to match with OSeMOSYS global technologies"""
PW2050_TECH_DICT = {'Hydro' : 'HYD', 
                    'Solar|CSP' : 'CSP',
                    'Solar|PV' : 'SPV',
                    'Wind|Onshore' : 'WON',
                    'Wind|Offshore' : 'WOF'}

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