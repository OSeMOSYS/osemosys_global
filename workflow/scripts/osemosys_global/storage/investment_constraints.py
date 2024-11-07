"""Function to set min and max capacity investment constraints."""

import pandas as pd
import itertools

from data import(
    get_years,
    get_years_build_rates)

from utils import apply_dtypes

def cap_investment_constraints_sto(storage_set, df_max_cap_invest_base,
                                   build_rates, no_investment_techs, 
                                   start_year, end_year, region_name):
    
    techs = storage_set['VALUE'].str[:3].unique()
    df_max_cap_invest_sto = df_max_cap_invest_base.copy()
    
    # Set max annual investments to 0 in case no expansion is allowed for specific storage technology.
    if any(x in no_investment_techs for x in techs) == True:
        no_invest_techs_sto = [y for y in techs if y in no_investment_techs]
        years = [get_years(start_year, end_year)]

        data = storage_set.loc[storage_set['VALUE'].str.startswith(tuple(no_invest_techs_sto))]
        
        data = pd.DataFrame(list(itertools.product([region_name], 
                                                   'PWR' + data['VALUE'].unique(), years, str(0))), 
                            columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"]
        )
        
        data =  data.explode('YEAR').reset_index(drop = True)
            
        df_max_cap_invest_sto = pd.concat([df_max_cap_invest_sto, data])
    
        # Filter no invest technologies out of user defined build rates.
        no_invest_techs_pwr = ['PWR' + tech for tech in no_investment_techs]
        build_rates = build_rates.loc[~build_rates['TECHNOLOGY'].str.startswith(tuple(no_invest_techs_pwr))]
    
    # Set storage technology and nodal specific max annual investments if defined.    
    if not build_rates.empty:
        data = build_rates.copy()
        
        data = data.apply(get_years_build_rates, axis = 1)
        
        data = data.explode(
            'YEAR').rename(columns = {'MAX_BUILD' : 'VALUE'})
        data['REGION'] = region_name
        
        df_max_cap_invest_sto = pd.concat([df_max_cap_invest_sto, data], join = 'inner')
   
    df_max_cap_invest_sto = apply_dtypes(df_max_cap_invest_sto, "TotalAnnualMaxCapacityInvestment")

    return df_max_cap_invest_sto