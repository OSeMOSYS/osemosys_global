"""Function to set min and max capacity investment constraints."""

import pandas as pd

from data import(
    get_years,
    get_years_build_rates
    )

from utils import apply_dtypes

def cap_investment_constraints_trn(df_iar_trn_final, df_max_cap_invest_base,
                                   build_rates, no_investment_techs, 
                                   start_year, end_year, region_name):
    
    # Set max annual investments to 0 in case no expansion is allowed.
    if 'TRN' in no_investment_techs:
        max_cap_invest_data = []
         
        max_cap_invest_techs = list(set(df_iar_trn_final.loc[
            df_iar_trn_final['TECHNOLOGY'].str.startswith('TRN')][
            'TECHNOLOGY'].tolist()))
        for tech in max_cap_invest_techs:
            for year in get_years(start_year, end_year): 
                max_cap_invest_data.append([region_name, tech, year, 0])
    
        df_max_cap_invest_trn = pd.DataFrame(max_cap_invest_data,
                                             columns = ['REGION', 'TECHNOLOGY', 
                                                        'YEAR', 'VALUE'])
        
        df_max_cap_invest_trn = pd.concat([df_max_cap_invest_base, 
                                           df_max_cap_invest_trn], join = 'inner')

    # Set pathway specific max annual investments if defined.
    elif not build_rates.empty:
        df_max_cap_invest_trn = build_rates.copy()
        
        df_max_cap_invest_trn = df_max_cap_invest_trn.apply(get_years_build_rates, 
                                                            axis = 1)
        
        df_max_cap_invest_trn = df_max_cap_invest_trn.explode(
            'YEAR').rename(columns = {'MAX_BUILD' : 'VALUE'})
        
        df_max_cap_invest_trn = df_max_cap_invest_trn.loc[df_max_cap_invest_trn["YEAR"].between(
            start_year, end_year)]

        df_max_cap_invest_trn['REGION'] = region_name
        
        df_max_cap_invest_trn = pd.concat([df_max_cap_invest_base, 
                                           df_max_cap_invest_trn], join = 'inner')
    
    else:
        df_max_cap_invest_trn = df_max_cap_invest_base.copy()
        
    df_max_cap_invest_trn = apply_dtypes(df_max_cap_invest_trn, 
                                         "TotalAnnualMaxCapacityInvestment")

    return df_max_cap_invest_trn