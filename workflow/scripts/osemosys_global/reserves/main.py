import pandas as pd
import os

from read import(
    import_technologies,
    import_fuels)

from reserve_margin import set_reserve_margin
from reserve_margin_tag_fuel import set_reserve_margin_fuels
from reserve_margin_tag_technology import set_reserve_margin_technologies

def main(tech_set_base: pd.DataFrame,
         fuel_set_base: pd.DataFrame):
    
    # CALL FUNCTIONS
    
    df_reserve_margin = set_reserve_margin(margins, start_year, 
                                           end_year, region_name)
    
    df_reserve_margin_tag_fuel = set_reserve_margin_fuels(fuel_set_base, start_year, 
                                                          end_year, region_name)
    
    df_reserve_margin_tag_tech = set_reserve_margin_technologies(margins_technologies,
                                                                 tech_set_base, 
                                                                 start_year, 
                                                                 end_year, 
                                                                 region_name)
        
    # OUTPUT CSV's
    
    
    df_reserve_margin.to_csv(os.path.join(output_data_dir, 'ReserveMargin.csv'), 
                             index = None)
    
    df_reserve_margin_tag_fuel.to_csv(os.path.join(output_data_dir, 
                                                   'ReserveMarginTagFuel.csv'), 
                             index = None)
    
    df_reserve_margin_tag_tech.to_csv(os.path.join(output_data_dir, 
                                                   'ReserveMarginTagTechnology.csv'), 
                             index = None) 

if __name__ == "__main__":
    
    if "snakemake" in globals():     
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        region_name = snakemake.params.region_name
        margins = snakemake.params.reserve_margin
        margins_technologies = snakemake.params.reserve_margin_technologies
        output_data_dir = snakemake.params.output_data_dir
        file_tech_set = f'{output_data_dir}/TECHNOLOGY.csv'
        file_fuel_set = f'{output_data_dir}/FUEL.csv'  
        
    # The below else statement defines variables if the 'transmission/main' script is to be run locally
    # outside the snakemake workflow. This is relevant for testing purposes only! User inputs when running 
    # the full workflow need to be defined in the config file. 

    else:         
        start_year = 2021
        end_year = 2050
        region_name = 'GLOBAL'
        margins = {'RM1': [10, 2025, 2029],
                   'RM2': [15, 2030, 2050]}
        
        margins_technologies = {
            'BIO' : 90,
            'CCG' : 90,
            'COA' : 90,
            'COG' : 50,
            'CSP' : 30,
            'GEO' : 90,
            'HYD' : 30,
            'OCG' : 90,
            'OIL' : 90,
            'OTH' : 90,
            'PET' : 90,
            'SPV' : 0,
            'URN' : 90,
            'WAS' : 90,
            'WAV' : 10,
            'WOF' : 10,
            'WON' : 10,
            'TRN' : 10,
            'SDS' : 69,
            'LDS' : 77}
        output_data_dir = 'results/data'
        file_tech_set = f'{output_data_dir}/TECHNOLOGY.csv'
        file_fuel_set = f'{output_data_dir}/FUEL.csv'

    # SET INPUT DATA
    
    tech_set_base = import_technologies(file_tech_set)
    fuel_set_base = import_fuels(file_fuel_set)  
    
    input_data = {
        "tech_set_base" : tech_set_base,
        "fuel_set_base" : fuel_set_base,
    }
    
    # CALL MAIN
    main(**input_data)