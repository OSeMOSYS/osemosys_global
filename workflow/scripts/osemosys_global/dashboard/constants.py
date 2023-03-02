"""Module for initialization constants"""

# Top level initializations
_PLOT_THEME = "seaborn"
_GEOGRAPHIC_SCOPE = "System"
_PLOT_TYPE = "Area"

# Map Initializaion Values
_MAP_THEME = "carto-positron"
_TRANSMISSION_LINE = "Hide"
_NODE_LOCATION = "Centroid"
_MAP_ELEMENTS_SIZE = 10

# Input Data Initialization 
_INPUT_DATA_CHOICE = "SpecifiedAnnualDemand"
_INPUT_DATA_TECH = "all"

# Result Data Initialization 
_RESULT_DATA_CHOICE = "ProductionByTechnologyAnnual"
_RESULT_DATA_TECH = "all"

# Transmission Data Initialization 
_TRANSMISSION_DATA_CHOICE = "TotalCapacityAnnual"

# Generator Abbreviations
TECHS_CONFIG = {
    "BIO":{
        "nicename": "Biomass",
        "color":"darkgreen",
    },
    "CCG":{
        "nicename": "Combined Cycle Natural Gas",
        "color":"lightcoral",
    },
    "COA":{
        "nicename": "Coal",
        "color":"black",
    },
    "COG":{
        "nicename": "Cogeneration",
        "color":"peru",
    },
    "CSP":{
        "nicename": "Concentrated Solar Power",
        "color":"wheat",
    },
    "ELC":{
        "nicename": "Electricity",
        "color":"gold",
    },
    "GAS":{
        "nicename": "Natural Gas",
        "color":"orange",
    },
    "GEO":{
        "nicename": "Geothermal",
        "color":"darkseagreen",
    },
    "HYD":{
        "nicename": "Hydroelectric",
        "color":"dodgerblue",
    },
    "INT":{
        "nicename":"International",
        "color":"darkgreen"
    },
    "OCG":{
        "nicename": "Open Cycle Natural Gas",
        "color":"firebrick",
    },
    "OIL":{
        "nicename": "Oil",
        "color":"lightgrey",
    },
    "OTH":{
        "nicename": "Other",
        "color":"teal",
    },
    "PET":{
        "nicename": "Petroleum",
        "color":"grey",
    },
    "SPV":{
        "nicename": "Solar PV",
        "color":"gold",
    },
    "URN":{
        "nicename": "Nuclear",
        "color":"mediumseagreen",
    },
    "WAV":{
        "nicename": "Wave",
        "color":"darkgreen",
    },
    "WAS":{
        "nicename": "Waste",
        "color":"navy",
    },
    "WOF":{
        "nicename": "Offshore Wind",
        "color":"violet",
    },
    "WON":{
        "nicename": "Onshore Wind",
        "color":"blueviolet",
    },
}

# Input Data
PARAM_CONFIG = {
    "SpecifiedAnnualDemand":{
        "nicename":"Annual Demand",
        "groupby":"FUEL",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"DEMAND (PJ)",
        "add_default":True,
        "default":0
    },
    "CapacityFactor":{
        "nicename":"Capacity Factor",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"mean",
        "xaxis":"TIMESLICE", 
        "ylabel":"",
        "add_default":False,
    },
    "CapitalCost":{
        "nicename":"Capital Costs",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",        
        "groupby_method":"mean",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M/GW)",
        "add_default":False,
    },
    "SpecifiedDemandProfile":{
        "nicename":"Demand Profile",
        "groupby":"FUEL",
        "groupby_method":"mean",
        "xaxis":"TIMESLICE", 
        "ylabel":"",
        "add_default":True,
        "default":0
    },
    "FixedCost":{
        "nicename":"Fixed Costs",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",        
        "groupby_method":"mean",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M/GW)",
        "add_default":False,
    },
    "ResidualCapacity":{
        "nicename":"Residual Capacity",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CAPACITY (GW)",
        "add_default":True,
        "default":0
    },
    "VariableCost":{
        "nicename":"Variable Costs",
        "filterby":"MIN",
        "groupby":"TECHNOLOGY",
        "groupby_method":"mean",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M/PJ)",
        "add_default":False,
    },
}

RESULT_CONFIG = {
    # "TotalDiscountedCost":{
    #     "nicename":"Discounted Cost",
    #     "groupby":None,
    #     "groupby_method":None,
    #     "xaxis":"YEAR", 
    #     "ylabel":"COST ($M)"
    # },
    "AnnualTechnologyEmission":{
        "nicename":"Emissions",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CO2 (MT)",
        "add_default":True,
        "default":0
    },
    "AnnualFixedOperatingCost":{
        "nicename":"Fixed Operating Cost",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M)",
        "add_default":True,
        "default":0
    },
    "NewCapacity":{
        "nicename":"New Capacity",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CAPACITY (GW)",
        "add_default":True,
        "default":0
    },
    "ProductionByTechnologyAnnual":{
        "nicename":"Production (Annual)",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"PRODUCTION (PJ)",
        "add_default":True,
        "default":0
    },
    "ProductionByTechnology":{
        "nicename":"Production (Time Slice)",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"TIMESLICE", 
        "ylabel":"PRODUCTION (PJ)",
        "add_default":True,
        "default":0
    },
    "TotalCapacityAnnual":{
        "nicename":"Total Capacity",
        "filterby":"PWR",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CAPACITY (GW)",
        "add_default":True,
        "default":0
    },
    "AnnualVariableOperatingCost":{
        "nicename":"Variable Operating Cost",
        "filterby":"MIN",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M)",
        "add_default":True,
        "default":0
    },
}

TRANSMISSION_CONFIG = {
    "NewCapacity":{
        "nicename":"New Capacity",
        "xaxis":"YEAR", 
        "ylabel":"CAPACITY (GW)",
        "add_default":True,
        "default":0
    },
    "ProductionByTechnologyAnnual":{
        "nicename":"Total Production (Annual)",
        "xaxis":"YEAR", 
        "ylabel":"PRODUCTION (PJ)",
        "add_default":True,
        "default":0
    },
    "ProductionByTechnology":{
        "nicename":"Total Production (Time Slice)",
        "xaxis":"TIMESLICE", 
        "ylabel":"PRODUCTION (PJ)",
        "add_default":True,
        "default":0
    },
    "TotalCapacityAnnual":{
        "nicename":"Total Capacity",
        "xaxis":"YEAR",
        "ylabel":"CAPACITY (GW)",
        "add_default":True,
        "default":0
    },
}