"""Function to calculate storage costs."""
import pandas as pd
import itertools

from data import(
    get_years
    )

def set_storage_capex_costs(storage_set, storage_param, 
                            start_year, end_year, 
                            region_name):

    capex_dict = {}
    duration_dict = {}
    
    # Set baseline technology capital costs and duration as defined in the config file.
    for tech, tech_params in storage_param.items():
        capex_dict[tech] = tech_params[0]
        duration_dict[tech] = tech_params[4]
    
    years = [get_years(start_year, end_year)]
    
    # CapitalCostStorage
    df_cap_cost_storage = pd.DataFrame(
        list(itertools.product([region_name], storage_set['VALUE'].unique(), years)), 
        columns=["REGION", "STORAGE", "YEAR"]
    )
    
    df_cap_cost_storage = df_cap_cost_storage.explode('YEAR').reset_index(drop = True)
    
    ''' Sets capital cost by taking the defined capital cost (m$/GW) divided by the storage 
    duration (=Storage Capacity (GWh)/Storage Power Rating (GW)) to get to GWh values followed 
    by the conversion to PJ (1 GWh = 0.0036 PJ).'''
    for tech, tech_params in storage_param.items():
        df_cap_cost_storage.loc[df_cap_cost_storage['STORAGE'].str.startswith(tech),
                   'VALUE'] = capex_dict[tech] / duration_dict[tech] / 0.0036
    
    return df_cap_cost_storage

def set_storage_operating_costs(storage_set, storage_param,
                                fom_base, var_base,
                                start_year, end_year, 
                                region_name):

    fom_dict = {}
    var_dict = {}
    duration_dict = {}
    
    # Set baseline costs and duration as defined in the config file.
    for tech, tech_params in storage_param.items():
        fom_dict[tech] = tech_params[1]
        var_dict[tech] = tech_params[2]
        duration_dict[tech] = tech_params[4]
    
    years = [get_years(start_year, end_year)]
    
    # FixedCost
    df_fom_storage = pd.DataFrame(
        list(itertools.product([region_name], storage_set['VALUE'].unique(), years)), 
        columns=["REGION", "TECHNOLOGY", "YEAR"]
    )
    
    df_fom_storage = df_fom_storage.explode('YEAR').reset_index(drop = True)
    
    
    ''' Sets fixed cost by taking the defined fixed cost (m$/GW/yr) divided by the storage 
    duration (=Storage Capacity (GWh)/Storage Power Rating (GW)) to mimic fixed costs for storage
    as capacities are set in GWh/PJ terms.'''
    for tech, tech_params in storage_param.items():
        df_fom_storage.loc[df_fom_storage['TECHNOLOGY'].str.startswith(tech),
                   'VALUE'] = fom_dict[tech] / duration_dict[tech]
        
    df_fom_storage['TECHNOLOGY'] = 'PWR' + df_fom_storage['TECHNOLOGY']
        
    df_fom_storage = pd.concat([fom_base, df_fom_storage])
    
    # VariableCost
    df_var_storage = pd.DataFrame(
        list(itertools.product([region_name], storage_set['VALUE'].unique(), str(2), years)), 
        columns=["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"]
    )
    
    df_var_storage = df_var_storage.explode('YEAR').reset_index(drop = True)
    
    ''' Sets variable cost by taking the defined variable cost ($/MWh) divided by the storage 
    duration (=Storage Capacity (GWh)/Storage Power Rating (GW)) to mimic variable costs for storage
    as capacities are set in GWh/PJ terms.'''
    for tech, tech_params in storage_param.items():
        df_var_storage.loc[df_var_storage['TECHNOLOGY'].str.startswith(tech),
                   'VALUE'] = round(var_dict[tech] / duration_dict[tech] / 3.6 , 4)
        
    df_var_storage['TECHNOLOGY'] = 'PWR' + df_var_storage['TECHNOLOGY']
        
    df_var_storage = pd.concat([var_base, df_var_storage])
    
    return df_fom_storage, df_var_storage
