"""Function to calculate residual capacity for powerplant technologies."""

import pandas as pd

from constants import (
    region_name,
    years
)

from utils import apply_dtypes

def cap_investment_constraints(df_iar_final, no_investment_techs):

    # Create totalAnnualMaxCapacityInvestment data 

    # Do not allow capacity investment for all PWRxxxxxxxx00 technolgoies 
    max_cap_invest_techs = list(set(
        df_iar_final.loc[df_iar_final['TECHNOLOGY'].str.endswith('00')]['TECHNOLOGY'].tolist()))
    max_cap_invest_data = []
    for tech in max_cap_invest_techs:
        for year in years: 
            max_cap_invest_data.append([region_name, tech, year, 0])

    # Do not allow investment for all xxxABCxxxxxxx technologies
    
    if not no_investment_techs:
        no_investment_techs = [] # Change from None type to empty list
    max_cap_invest_techs = list(set(df_iar_final.loc[
        df_iar_final['TECHNOLOGY'].str[3:6].isin(no_investment_techs)][
        'TECHNOLOGY'].tolist()))
    for tech in max_cap_invest_techs:
        for year in years: 
            max_cap_invest_data.append([region_name, tech, year, 0])
    
    # Save totalAnnualMaxCapacityInvestment
    df_max_cap_invest = pd.DataFrame(max_cap_invest_data,
                                    columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )
    df_max_cap_invest = apply_dtypes(df_max_cap_invest, "TotalAnnualMaxCapacityInvestment")

    df_min_cap_invest = pd.DataFrame(columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )
    df_min_cap_invest = apply_dtypes(df_min_cap_invest, "TotalAnnualMinCapacityInvestment")

    return df_max_cap_invest, df_min_cap_invest