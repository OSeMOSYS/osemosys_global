"""Function to calculate residual capacity for powerplant technologies."""

import pandas as pd
import numpy as np
from utils import apply_dtypes

def costs_pwr(df_weo_data, costs_dict):

    # ### Costs: Capital, fixed, and variable

    df_costs = pd.melt(df_weo_data, 
                       id_vars = ['technology', 'weo_region', 'parameter'], 
                       value_vars = ['2019', '2030', '2040'], 
                       var_name = 'YEAR')
    df_costs['parameter'] = df_costs['parameter'].str.split('\r\n').str[0]
    df_costs['value'] = df_costs['value'].replace({'n.a.':0})
    df_costs['value'] = df_costs['value'].astype(float) 
    df_costs = df_costs.pivot_table(index = ['technology', 'parameter', 'YEAR'], 
                                    columns = 'weo_region', 
                                    values = 'value').reset_index()
    df_costs['AS_average'] = (df_costs['China'] + 
                                df_costs['India'] + 
                                df_costs['Japan'] + 
                                df_costs['Middle East']).div(4)
    df_costs['SEA_average'] = (df_costs['Indonesia'] + 
                                df_costs['Vietnam']).div(2)
    df_costs['NA_average'] = (df_costs['United States'])
    df_costs['SA_average'] = (df_costs['Brazil'])
    df_costs['Global_average'] = (df_costs['Africa'] +
                                  df_costs['Brazil'] + 
                                  df_costs['Europe'] +
                                  df_costs['China'] + 
                                  df_costs['India'] + 
                                  df_costs['Japan'] + 
                                  df_costs['Middle East'] +
                                  df_costs['Russia'] +
                                  df_costs['United States']).div(9)
    df_costs = pd.melt(df_costs, 
                       id_vars = ['technology', 'parameter', 'YEAR'], 
                       value_vars = [x 
                                     for x 
                                     in df_costs.columns 
                                     if x not in ['technology', 'parameter', 'YEAR']
                                    ]
                      )
    df_costs['YEAR'] = df_costs['YEAR'].astype(int)

    df_costs = df_costs.loc[df_costs['technology'].isin(costs_dict.keys())]
    df_costs['technology_code'] = df_costs['technology'].replace(costs_dict)
    
    return df_costs

def costs_end(df_weo_regions, df_costs, df_oar_final, df_trans_capex, 
              df_trans_fix):    
    
    weo_regions_dict = dict([(k, v) 
                             for k, v 
                             in zip(df_weo_regions['technology_code'], 
                                    df_weo_regions['weo_region']
                                   )
                            ]
                           )

    # Create formatted CSVs
    for each_cost in ['Capital', 'O&M']:
        df_costs_temp = df_costs.copy().loc[df_costs.copy()['parameter'].str.contains(each_cost)]
        
        df_costs_temp.drop(['technology', 'parameter'], 
                           axis = 1, 
                           inplace = True)
        df_costs_final = df_oar_final[['REGION',
                                       'TECHNOLOGY',
                                       'YEAR'
                                      ]]
        df_costs_final.loc[:,'YEAR'] = df_costs_final['YEAR'].astype(int)
        
        df_costs_final = df_costs_final.drop_duplicates()
        df_costs_final = (df_costs_final
                          .loc[(df_costs_final['TECHNOLOGY']
                                .str.startswith('PWR')
                               ) & 
                               (~df_costs_final['TECHNOLOGY']
                                .str.contains('TRN')
                               )
                              ]
                         )
        df_costs_final['technology_code'] = df_costs_final['TECHNOLOGY'].str[3:6]
        df_costs_final['weo_region'] = df_costs_final['TECHNOLOGY'].str[6:9]
        df_costs_final['weo_region'] = (df_costs_final['weo_region']
                                             .replace(weo_regions_dict))

        df_costs_final = pd.merge(df_costs_final, 
                                  df_costs_temp, 
                                  on = ['technology_code', 'weo_region', 'YEAR'], 
                                  how = 'left'
                                 )
        df_costs_final.drop(['technology_code', 'weo_region'], 
                            axis = 1, 
                            inplace = True)
        df_costs_final = df_costs_final.infer_objects().fillna(-9)
        df_costs_final = pd.pivot_table(df_costs_final, 
                                        index = ['REGION', 'YEAR'], 
                                        columns = 'TECHNOLOGY', 
                                        values = 'value').reset_index()
        df_costs_final = df_costs_final.replace([-9],[np.nan])

        for col in df_costs_final.columns:
            if df_costs_final[col].dtype != 'object':
                df_costs_final[col] = df_costs_final[col].interpolate(method = 'linear', 
                                                            limit_direction='forward').round(2)
                
                df_costs_final[col] = df_costs_final[col].interpolate(method = 'linear', 
                                                            limit_direction='backward').round(2)

        df_costs_final = pd.melt(df_costs_final, 
                                 id_vars = ['REGION', 'YEAR'], 
                                 value_vars = [x for x in df_costs_final.columns
                                               if x not in ['REGION', 'YEAR']
                                              ],
                                 var_name = 'TECHNOLOGY',
                                 value_name = 'VALUE'
                                )
        df_costs_final = df_costs_final[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]
        df_costs_final = df_costs_final[~df_costs_final['VALUE'].isnull()]

        if each_cost in ['Capital']:
            df_costs_final_capital = df_costs_final.merge(df_trans_capex, how='outer')
            df_costs_final_capital = apply_dtypes(df_costs_final_capital, "CapitalCost")

            
        if each_cost in ['O&M']:
            df_costs_final_fixed = df_costs_final.merge(df_trans_fix, how='outer')
            df_costs_final_fixed = apply_dtypes(df_costs_final_fixed, "FixedCost")
       
    return df_costs_final_capital, df_costs_final_fixed