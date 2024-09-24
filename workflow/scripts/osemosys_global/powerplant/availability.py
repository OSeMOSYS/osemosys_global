"""Function to set availability factors."""

import pandas as pd
import itertools

from constants import(
    years,
    region_name
    )

def availability_factor(availability, tech_set_base):
    
    af_dict = dict(zip(list(availability['technology']),
                       list(availability['value'])))
    
    tech_list = [x for x in tech_set_base['VALUE']
                 if x.startswith('PWR')]
    df_af_final = pd.DataFrame(list(itertools.product(tech_list,
                                                      years)
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