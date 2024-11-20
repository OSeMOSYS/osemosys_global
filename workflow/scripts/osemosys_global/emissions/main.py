import pandas as pd
import os

from read import(
    import_emission_factors,
    import_iar_base,
    import_oar_base,
    get_ember_emissions
)

from constants import(
    _TECH_TO_FUEL,
    _EMISSION,
    CCS_EFF
    )

from emission_activity_ratio import get_ear

from sets import set_unique_emissions

from emission_penalty import get_emission_penalty

from emission_limit import add_emission_limits

def main(
        ember: pd.DataFrame,
        emission_factors: pd.DataFrame,
        emission_penalty: list,
        emission_limit: list, 
        iar_base: pd.DataFrame,
        oar_base: pd.DataFrame,        
):
    
    # CALL FUNCTIONS
    
    # Set EmissionActivityRatio.
    df_emission_activity_ratio = get_ear(_EMISSION, emission_factors, CCS_EFF,
                                         iar_base, oar_base, _TECH_TO_FUEL)
    
    # Set EMISSION set.
    df_emissions_set = set_unique_emissions(df_emission_activity_ratio)

    # Set user defined emission penalties.
    df_emission_penalty = get_emission_penalty(df_emissions_set, emission_penalty, 
                                               start_year, end_year, region_name)

    
    # Set user defined emission limits.
    df_annual_emission_limit = add_emission_limits(df_emissions_set, emission_limit, 
                                                   ember, start_year, end_year, 
                                                   region_name)

    # OUTPUT CSV's
    df_emission_activity_ratio.to_csv(os.path.join(output_data_dir, "EmissionActivityRatio.csv"), 
                           index = None)
    
    df_emissions_set.to_csv(os.path.join(output_data_dir, "EMISSION.csv"), 
                         index = None)
    
    df_emission_penalty.to_csv(os.path.join(output_data_dir, "EmissionsPenalty.csv"), 
                         index = None)
    
    df_annual_emission_limit.to_csv(os.path.join(output_data_dir, "AnnualEmissionLimit.csv"), 
                         index = None)


if __name__ == "__main__":
    
    if "snakemake" in globals():
        file_ember = snakemake.input.ember
        file_emission_factors = snakemake.input.emissions_factors         
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        region_name = snakemake.params.region_name
        output_data_dir = snakemake.params.output_data_dir
        input_data_dir = snakemake.params.input_data_dir
        file_iar_base = snakemake.input.iar
        file_oar_base = snakemake.input.oar
        emission_penalty = snakemake.params.emission_penalty
        emission_limit = snakemake.params.emission_limit
            
    # The below else statement defines variables if the 'powerplant/main' script is to be run locally
    # outside the snakemake workflow. This is relevant for testing purposes only! User inputs when running 
    # the full workflow need to be defined in the config file. 
    else:
        file_ember = 'resources/data/ember_yearly_electricity_data.csv'
        file_emission_factors = 'resources/data/emission_factors.csv'        
        start_year = 2021
        end_year = 2050
        region_name = 'GLOBAL'
        output_data_dir = 'results/data'
        input_data_dir = 'resources/data'
        file_iar_base = f'{output_data_dir}/InputActivityRatio.csv'
        file_oar_base = f'{output_data_dir}/OutputActivityRatio.csv'
        emission_penalty = [["CO2", "IND", 2020, 2050, 2.1]]
        emission_limit = [["CO2", "IND", "POINT", 2048, 0],
                          ["CO2", "IND", "LINEAR", 2040, 1],
                          ["CO2", "IND", "LINEAR", 2028, 400],
                          ["CO2", "IND", "POINT", 2030, 300]]

    # SET INPUT DATA
    ember = get_ember_emissions(file_ember)
    emission_factors = import_emission_factors(file_emission_factors)
    
    iar_base = import_iar_base(file_iar_base)
    oar_base = import_oar_base(file_oar_base)
    
    input_data = {
    'ember' : ember,
    'emission_factors' : emission_factors,
    'emission_penalty' : emission_penalty,
    'emission_limit' : emission_limit,
    'iar_base' : iar_base,
    'oar_base' : oar_base,
    }
    
    # CALL MAIN
    main(**input_data)