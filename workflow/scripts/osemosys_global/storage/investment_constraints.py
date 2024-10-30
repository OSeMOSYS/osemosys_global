"""Function to set min and max capacity investment constraints."""

import pandas as pd
import itertools

from data import get_years

from utils import apply_dtypes

def cap_investment_constraints_sto(storage_set, max_cap_invest_base,
                                   no_investment_techs, start_year, 
                                   end_year, region_name):
    
    techs = storage_set['VALUE'].str[:3].unique()
    if any(x in no_investment_techs for x in techs) == True:
        no_invest_techs_sto = [y for y in techs if y in no_investment_techs]
        years = [get_years(start_year, end_year)]

        data = storage_set.loc[storage_set['VALUE'].str.startswith(tuple(no_invest_techs_sto))]
        
        df_max_cap_invest_sto = pd.DataFrame(
            list(itertools.product([region_name], 'PWR' + data['VALUE'].unique(), years, str(0))), 
            columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"]
        )
        
        df_max_cap_invest_sto =  df_max_cap_invest_sto.explode('YEAR').reset_index(drop = True)
        
        df_max_cap_invest_sto = apply_dtypes(df_max_cap_invest_sto, "TotalAnnualMaxCapacityInvestment")
            
        df_max_cap_invest_sto = pd.concat([max_cap_invest_base, df_max_cap_invest_sto])
            
    else:
        df_max_cap_invest_sto = max_cap_invest_base.copy()

    return df_max_cap_invest_sto