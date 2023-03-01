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
_INPUT_DATA_TECH = "BIO"

# Result Data Initialization 
_RESULT_DATA_CHOICE = "AnnualTechnologyEmission"
_RESULT_DATA_TECH = "COA"

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
    "GEO":{
        "nicename": "Geothermal",
        "color":"darkseagreen",
    },
    "HYD":{
        "nicename": "Hydroelectric",
        "color":"dodgerblue",
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
        "ylabel":"DEMAND (PJ)"
    },
    "CapacityFactor":{
        "nicename":"Capacity Factor",
        "groupby":"TECHNOLOGY",
        "groupby_method":"mean",
        "xaxis":"TIMESLICE", 
        "ylabel":""
    },
    "CapitalCost":{
        "nicename":"Capital Costs",
        "groupby":"TECHNOLOGY",        
        "groupby_method":"mean",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M/GW)"
    },
    "SpecifiedDemandProfile":{
        "nicename":"Demand Profile",
        "groupby":"FUEL",
        "groupby_method":"mean",
        "xaxis":"TIMESLICE", 
        "ylabel":""
    },
    "FixedCost":{
        "nicename":"Fixed Costs",
        "groupby":"TECHNOLOGY",        
        "groupby_method":"mean",
        "xaxis":"YEAR", 
        "ylabel":"COOST ($M/GW)"
    },
    "ResidualCapacity":{
        "nicename":"Residual Capacity",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CAPACITY (GW)"
    },
    "VariableCost":{
        "nicename":"Variable Costs",
        "groupby":"TECHNOLOGY",
        "groupby_method":"mean",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M/PJ)"
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
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CO2 (MT)"
    },
    "AnnualFixedOperatingCost":{
        "nicename":"Fixed Operating Cost",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M)"
    },
    "NewCapacity":{
        "nicename":"New Capacity",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CAPACITY (GW)"
    },
    "ProductionByTechnologyAnnual":{
        "nicename":"Production (Annual)",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"PRODUCTION (PJ)"
    },
    "ProductionByTechnology":{
        "nicename":"Production (Time Slice)",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"TIMESLICE", 
        "ylabel":"PRODUCTION (PJ)"
    },
    "TotalCapacityAnnual":{
        "nicename":"Total Capacity",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "xaxis":"YEAR", 
        "ylabel":"CAPACITY (GW)"
    },
    "AnnualVariableOperatingCost":{
        "nicename":"Variable Operating Cost",
        "groupby":"TECHNOLOGY",
        "groupby_method":"sum",
        "parsing_method":"MIN",
        "xaxis":"YEAR", 
        "ylabel":"COST ($M)"
    },
}