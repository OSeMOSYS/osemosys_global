"""Rules for validating historical results"""

###
# constants
###

# used for target rules
CAPACITY_VALIDATION = ["ember", "irena", "eia"]
GENERATION_VALIDATION = ["ember", "irena", "eia"]
EMISSION_VALIDATION = ["ember", "climatewatch"]
EMISSION_INTENSITY_VALIDATION = ["ember"]

###
# capacity
###

def capacity_validation_data(wildcards):
    if wildcards.datasource == "ember":
        return "resources/data/ember_yearly_electricity_data.csv"
    elif wildcards.datasource == "irena":
        return "resources/data/validation/irena_capacity.csv"
    elif wildcards.datasource == "eia":
        return "resources/data/validation/eia_capacity.json"

def capacity_options(wildcards):
    if wildcards.datasource == "irena":
        return {"iso_codes": "resources/data/validation/iso.csv"}
    else:
        return []

rule validate_capacity:
    message: "Validating capacity against {wildcards.datasource}"
    input:
        unpack(capacity_options),
        validation_data = capacity_validation_data,
        og_result = "results/{scenario}/results/TotalCapacityAnnual.csv",
    params:
        result_dir="results/{scenario}/results",
        variable="capacity"
    output:
        expand("results/{{scenario}}/validation/{country}/capacity/{{datasource}}.png", country=COUNTRIES)
    script:
        "../scripts/osemosys_global/validation/main.py"

###
# generation
###

def generation_validation_data(wildcards):
    if wildcards.datasource == "ember":
        return "resources/data/ember_yearly_electricity_data.csv"
    elif wildcards.datasource == "irena":
        return "resources/data/validation/irena_generation.csv"
    elif wildcards.datasource == "eia":
        return "resources/data/validation/eia_generation.json"

def generation_options(wildcards):
    if wildcards.datasource == "irena":
        return {"iso_codes": "resources/data/validation/iso.csv"}
    else:
        return []

rule validate_generation:
    message: "Validating generation against {wildcards.datasource}"
    input:
        unpack(generation_options),
        validation_data = generation_validation_data,
        og_result = "results/{scenario}/results/ProductionByTechnologyAnnual.csv",
    params:
        result_dir="results/{scenario}/results",
        variable="generation"
    output:
        expand("results/{{scenario}}/validation/{country}/generation/{{datasource}}.png", country=COUNTRIES)
    script:
        "../scripts/osemosys_global/validation/main.py"

###
# emissions
###

def emission_validation_data(wildcards):
    if wildcards.datasource == "ember":
        return "resources/data/ember_yearly_electricity_data.csv"
    elif wildcards.datasource == "climatewatch":
        return "resources/data/validation/climatewatch.csv"

rule validate_emissions:
    message: "Validating emissions against {wildcards.datasource}"
    input:
        validation_data = emission_validation_data,
        og_result = "results/{scenario}/results/AnnualEmissions.csv"
    params:
        result_dir="results/{scenario}/results",
        variable="emissions"
    output:
        expand("results/{{scenario}}/validation/{country}/emissions/{{datasource}}.png", country=COUNTRIES)
    script:
        "../scripts/osemosys_global/validation/main.py"

###
# emission intensity
###

def emission_intensity_validation_data(wildcards):
    if wildcards.datasource == "ember":
        return "resources/data/ember_yearly_electricity_data.csv"

rule validate_emission_intensity:
    message: "Validating emission intensity against {wildcards.datasource}"
    input:
        validation_data = emission_intensity_validation_data,
        og_result = "results/{scenario}/results/AnnualEmissionIntensity.csv"
    params:
        result_dir="results/{scenario}/results",
        variable="emission_intensity"
    output:
        expand("results/{{scenario}}/validation/{country}/emission_intensity/{{datasource}}.png", country=COUNTRIES)
    script:
        "../scripts/osemosys_global/validation/main.py"