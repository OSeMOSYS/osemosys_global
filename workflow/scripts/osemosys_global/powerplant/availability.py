"""Function to set availability factors."""

import pandas as pd
import itertools

from data import get_years

def availability_factor(availability, tech_set_base,
                        start_year, end_year, region_name):
    
    af_dict = dict(zip(list(availability['technology']),
                       list(availability['value'])))
    
    tech_list = [x for x in tech_set_base['VALUE']
                 if x.startswith('PWR')]
    df_af_final = pd.DataFrame(list(itertools.product(tech_list,
                                                      get_years(start_year, end_year))
                                    ),
                               columns = ['TECHNOLOGY', 'YEAR']
                               )
    df_af_final['TECH'] = df_af_final['TECHNOLOGY'].str[3:6]
    df_af_final['VALUE'] = df_af_final['TECH'].map(af_dict)
    df_af_final.dropna(inplace=True)
    df_af_final['REGION'] = region_name
    df_af_final = df_af_final[['REGION',
                               'TECHNOLOGY',
                               'YEAR',
                               'VALUE']]
    
    return df_af_final

def availability_factor_custom(df_af, custom_af_factors):
    
    if custom_af_factors:
        
        # Loop through country level targets and apply
        for param in custom_af_factors:
            region = param[0]
            technology = param[1]
            first_year = param[2]
            last_year = param[3]
            value = param[4]
            
            if len(region) == 3:
                df_af.loc[(df_af["YEAR"].between(first_year, last_year)) & 
                          (df_af["TECHNOLOGY"].str[6:9] == region) & 
                          (df_af["TECHNOLOGY"].str[3:6] == technology) 
                          ,"VALUE"] = value / 100
        
        # Do the same for nodal level targets to overwrite country level targets
        for param in custom_af_factors:
            region = param[0]
            technology = param[1]
            first_year = param[2]
            last_year = param[3]
            value = param[4]
            
            if len(region) == 5:
                df_af.loc[(df_af["YEAR"].between(first_year, last_year)) & 
                          (df_af["TECHNOLOGY"].str[6:11] == region) & 
                          (df_af["TECHNOLOGY"].str[3:6] == technology) 
                          ,"VALUE"] = value / 100
    
    return df_af