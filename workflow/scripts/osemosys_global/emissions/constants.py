"""Constants for the emissions module"""

# Constant for tech to fuel emission type mapping
_TECH_TO_FUEL = {
    #'BIO':'Bagasse',
    "WAS": "Municipal Solid Waste",
    "COA": "Lignite Coal",
    "COG": "Lignite Coal",
    "OCG": "Natural Gas",
    "CCG": "Natural Gas",
    "GAS": "Natural Gas",
    "PET": "Crude Oil",
    "OIL": "Crude Oil",
    "OTH": "Natural Gas",
    "CCS": "Lignite Coal",
}

# Emission name
_EMISSION = "CO2"

# COA CCS CO2 capture efficiency
CCS_EFF = 90