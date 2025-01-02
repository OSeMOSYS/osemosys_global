"""Constants for the storage module"""

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

# Group GESDB technologies into OG technologies.
GESDB_TECH_MAP = {'Compressed air energy storage' : 'SDS', 
                  'Electro-chemical capacitor' : 'SDS', 
                  'Flow battery' : 'SDS', 
                  'Flywheel' : 'SDS', 
                  'Heat thermal storage' : 'SDS', 
                  'Hydrogen storage' : 'SDS',
                  'Latent heat' : 'SDS', 
                  'Lead-acid battery' : 'SDS', 
                  'Lithium-ion battery' : 'SDS', 
                  'Nickel-based battery' : 'SDS', 
                  'Pumped hydro storage' : 'LDS', 
                  'Sensible heat' : 'SDS', 
                  'Sodium-based battery' : 'SDS', 
                  'Zinc-based battery' : 'SDS',
                  'null' : 'SDS'}

'''Set duration type to 'Default' for standardized values following the user defined 
'storage_parameters' from the config file or 'Historical' for node and aggregate technology 
specific duration values (i.e. the ratio between storage capacity (kWh) and power ratings (kW)) 
as derived from the GESDB datasets. If 'Historical' is chosen the entries that do not have 
assigned storage capacity values (kWh) in the original dataset will still use the 
'storage_parameters' defined duration.'''
DURATION_TYPE = 'Historical'

# Build year for plants without commissioning year attached in dataset
BUILD_YEAR = 2030
# Retirement year for plants that are already past its assumed operational lifetime
RETIREMENT_YEAR = 2030

INACTIVE_VARS = ['Decomissioned', 'Decommissioned', 'De-Commissioned', ]

ACTIVE_VARS = ['Operational', 'Standby/backup', 'Offline/Under Repair', 'Offline/under repair',]

NEW_VARS = ['Analysis', 'Analysis Only', 'Annocunced/never built', 'Announced', 
            'Announced/never built', 'Announced/Never Built', 'Constructed', 
            'Contracted', 'null', 'Procurement', 'Project Development', 
            'Under Construction', 'Under construction']
