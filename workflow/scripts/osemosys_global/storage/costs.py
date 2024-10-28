"""Function to calculate storage costs."""
import pandas as pd
import itertools

from data import(
    get_years
    )

def set_storage_costs(storage_set, storage_param, 
                      start_year, end_year, 
                      region_name):

    capex_dict = {}
    
    # Set baseline technology capital costs as defined in the config file.
    for tech, tech_params in storage_param.items():
        capex_dict[tech] = tech_params[0]
    
    years = [get_years(start_year, end_year)]
    
    # CapitalCostStorage
    df_cap_cost_storage = pd.DataFrame(
        list(itertools.product([region_name], storage_set['VALUE'].unique(), years)), 
        columns=["REGION", "STORAGE", "YEAR"]
    )
    
    df_cap_cost_storage = df_cap_cost_storage.explode('YEAR').reset_index(drop = True)
    
    for tech, tech_params in storage_param.items():
        df_cap_cost_storage.loc[df_cap_cost_storage['STORAGE'].str.startswith(tech),
                   'VALUE'] = capex_dict[tech]
    
    return df_cap_cost_storage