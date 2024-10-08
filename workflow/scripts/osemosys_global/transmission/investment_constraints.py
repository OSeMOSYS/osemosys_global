"""Function to set min and max capacity investment constraints."""

import pandas as pd

from data import get_years

from utils import apply_dtypes

def cap_investment_constraints_trn(df_iar_trn_final, df_max_cap_invest_base,
                                   no_investment_techs, start_year, 
                                   end_year, region_name):

    if 'TRN' in no_investment_techs:
        max_cap_invest_data = []
         
        max_cap_invest_techs = list(set(df_iar_trn_final.loc[
            df_iar_trn_final['TECHNOLOGY'].str.startswith('TRN')][
            'TECHNOLOGY'].tolist()))
        for tech in max_cap_invest_techs:
            for year in get_years(start_year, end_year): 
                max_cap_invest_data.append([region_name, tech, year, 0])
    
        # Save totalAnnualMaxCapacityInvestment
        df_max_cap_invest_trn = pd.DataFrame(max_cap_invest_data,
                                        columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                        )
        df_max_cap_invest_trn = apply_dtypes(df_max_cap_invest_trn, "TotalAnnualMaxCapacityInvestment")
        
        df_max_cap_invest_trn = pd.concat([df_max_cap_invest_base, df_max_cap_invest_trn])
        
    else:
        df_max_cap_invest_trn = df_max_cap_invest_base

    return df_max_cap_invest_trn