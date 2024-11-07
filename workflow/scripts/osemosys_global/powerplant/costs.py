"""Function to calculate powerplant technology costs."""

import pandas as pd
import numpy as np
import itertools

from data import get_years

from constants import COSTS_DICT

from utils import apply_dtypes


def costs_pwr(df_weo_data):

    # ### Costs: Capital, fixed

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

    df_costs = df_costs.loc[df_costs['technology'].isin(COSTS_DICT.keys())]
    df_costs['technology_code'] = df_costs['technology'].replace(COSTS_DICT)
    
    return df_costs

def costs_end(df_weo_regions, df_costs, df_oar_final): 
    
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
            df_costs_final_capital = apply_dtypes(df_costs_final, "CapitalCost")

            
        if each_cost in ['O&M']:
            df_costs_final_fixed = apply_dtypes(df_costs_final, "FixedCost")
       
    return df_costs_final_capital, df_costs_final_fixed

def calculate_cmo_forecast_prices(df_cmo_data, 
                                  cmo_data_year,
                                  biomass_var_costs,
                                  nuclear_var_costs,
                                  waste_var_costs,
                                  int_cost_factor):
    
    '''Data for variable costs of fuels taken from World Bank Commodity 
    Market Outlooks: https://www.worldbank.org/en/research/commodity-markets'''

    data = df_cmo_data[['Commodity', cmo_data_year]].replace(
        {'Coal, Australia' : 'COA',
         'Crude oil, Brent' : 'OIL',
         'Natural gas, Europe' : 'KFNGAS_EU',
         'Natural gas, U.S.' : 'KFNGAS_US',
         'Liquefied natural gas, Japan' : 'KFNGAS_JP',
         }, regex = True).rename(columns = {'Commodity' : 'FUEL',
                                             cmo_data_year : 'VALUE'})
    
    '''Add average GAS price based on the 'KFNGAS_EU', KFNGAS_US' and
    KFNGAS_JP' entries.'''
    gas_avg = data.loc[data['FUEL'].str.contains('GAS')].loc[:, 'VALUE'].mean()
    data.loc[len(data.index)] = ['GAS', gas_avg]
    data = data.loc[~data['FUEL'].isin(['KFNGAS_EU', 'KFNGAS_US', 'KFNGAS_JP'])]
    
    '''Convert to $/PJ
    # Values taken from kylesconverter.com:
    #  1 mt coal contains 29.31 GJ. We want $mill/PJ so divide by 29.31
    #  1 mmbtu = 0.000001055056 PJ. We want $mill/PJ so divide by 0.000001055056 
    and multiply with 1000000 (= / 1.055056).
    #  1 bbl = 0.00000612 PJ (Barrels of Oil) Original was $/bbl so divide by 
    0.00000612 and divide by 1000000 ( = / 6.12)
    
    # ORIGINAL UNITS:
    #       MINCOA MINOIL MINGAS  KFNGAS_US KFNGAS_JP
    #Unit   $/mt   $/bbl  $/mmbtu $/mmbtu   $/mmbtu'''
    data.loc[data['FUEL'] == 'COA', 'VALUE'] = data['VALUE']  / 29.31
    data.loc[data['FUEL'] == 'OIL', 'VALUE'] = data['VALUE']  / 6.12
    data.loc[data['FUEL'] == 'GAS', 'VALUE'] = data['VALUE']  / 1.055056
    
    '''Add in other fuels that are the same as those above. Cogen is powered by gas.
    Other petroleum products are similar to oil. Petroleum products are similar to oil.'''
    
    data = pd.concat([data, data.loc[data['FUEL'] == 'GAS'].replace('GAS', 'COG')])
    data = pd.concat([data, data.loc[data['FUEL'] == 'OIL'].replace('OIL', 'OTH')])
    data = pd.concat([data, data.loc[data['FUEL'] == 'OIL'].replace('OIL', 'PET')])

    # Add in default costs.
    data.loc[len(data.index)] = ['BIO', biomass_var_costs] 
    data.loc[len(data.index)] = ['URN', nuclear_var_costs]
    data.loc[len(data.index)] = ['WAS', waste_var_costs]
    
    data['COUNTRY'] = 'INT'
    data['UNIT'] = 'm$/PJ'
    data['ENERGY_CONTENT'] = 1
    
    # Multiply fuel prices with international cost factor.
    data['VALUE'] = data['VALUE'] * int_cost_factor
    
    # Set costs for all years equal to custom fuel prices csv.
    for year in [2020, 2025, 2030, 2040, 2050]:
        data[str(year)] = data['VALUE']
        
    data = data.drop(columns = {'VALUE'}).reset_index(drop = True)

    return data

def costs_var_pwr(df_fuel_prices, tech_set,
                  start_year, end_year, region_name,
                  df_cmo_data, cmo_data_year,
                  biomass_var_costs, nuclear_var_costs,
                  waste_var_costs, int_cost_factor):
    
    # Pulls CMO data and custom inputs from constants.
    cmo_inputs = calculate_cmo_forecast_prices(df_cmo_data, 
                                      cmo_data_year,
                                      biomass_var_costs,
                                      nuclear_var_costs,
                                      waste_var_costs,
                                      int_cost_factor)
    
    # Checks if INT prices are user defined and if not adds entries.        
    for fuel in cmo_inputs['FUEL'].unique():
        if df_fuel_prices.loc[(df_fuel_prices['FUEL'] == fuel) & 
                              (df_fuel_prices['COUNTRY'] == 'INT')].empty:
            data = cmo_inputs.loc[cmo_inputs['FUEL'] == fuel]
            df_fuel_prices = pd.concat([df_fuel_prices, data]
                                       ).reset_index(drop = True)
        
    # Read in Technologies
    years =  get_years(start_year, end_year)
    # ### Filter technologies to keep only fuel production technologies (MIN)

    df_techs = tech_set[tech_set.VALUE.str.contains('|'.join(['MIN','RNW']))]
    techs = df_techs['VALUE']

    # New calculation of fuel prices
    df_fuel_prices.drop(['UNIT'],
                        axis=1,
                        inplace=True)
    df_fuel_prices = pd.melt(df_fuel_prices,
                             id_vars=['FUEL',
                                      'COUNTRY',
                                      'ENERGY_CONTENT'],
                             value_vars=[x for x in df_fuel_prices.columns
                                         if x not in ['FUEL',
                                                      'COUNTRY',
                                                      'ENERGY_CONTENT']],
                             var_name='YEAR',
                             value_name='VALUE')
    df_fuel_prices['VALUE'] = (df_fuel_prices['VALUE'] / 
                               df_fuel_prices['ENERGY_CONTENT'])
    #df_fuel_prices['TECHNOLOGY'] = ('MIN' +
    #                                df_fuel_prices['FUEL'] +
    #                                df_fuel_prices['COUNTRY'])
    df_fuel_prices.loc[df_fuel_prices['FUEL'].isin(['BIO','WAS']),
                       'TECHNOLOGY'] = ('RNW' +
                                        df_fuel_prices['FUEL'] +
                                        df_fuel_prices['COUNTRY'])      
    df_fuel_prices.loc[~(df_fuel_prices['FUEL'].isin(['BIO','WAS'])),
                       'TECHNOLOGY'] = ('MIN' +
                                        df_fuel_prices['FUEL'] +
                                        df_fuel_prices['COUNTRY'])
    df_fuel_prices['VALUE'] = df_fuel_prices['VALUE'].round(2)

    # Dataframe with INT fuel prices
    df_int_fuel_prices = df_fuel_prices.loc[df_fuel_prices['COUNTRY'].isin(['INT'])]
    df_int_fuel_prices = df_int_fuel_prices[['FUEL',
                                             'YEAR',
                                             'VALUE']]
    df_int_fuel_prices['YEAR'] = df_int_fuel_prices['YEAR'].astype(int)

    # List of countries with country-specific values
    country_list = [x for x in df_fuel_prices['COUNTRY'].unique()
                    if x not in ['INT']]

    # Dataframe with country-specific and INT values
    df_fuel_prices = df_fuel_prices[['TECHNOLOGY',
                                     'YEAR',
                                     'VALUE']]
    df_fuel_prices['YEAR'] = df_fuel_prices['YEAR'].astype(int)
    df_fuel_prices.rename(columns={'TECHNOLOGY': 'TECH_COUNTRY'},
                          inplace=True)

    # Scaffolding for final dataframe
    df_fuel_prices_final = pd.DataFrame(list(itertools.product(list(techs),
                                                               [1, 2],
                                                               years)),
                                        columns = ['TECHNOLOGY',
                                                   'MODE_OF_OPERATION',
                                                   'YEAR']
                                        )
    df_fuel_prices_final['YEAR'] = df_fuel_prices_final['YEAR'].astype(int)
    df_fuel_prices_final['FUEL'] = df_fuel_prices_final['TECHNOLOGY'].str[3:6]
    df_fuel_prices_final['TECH_COUNTRY'] = df_fuel_prices_final['TECHNOLOGY'].str[:9]

    # Values for countries WITHOUT country-specific data set to INT values
    df_fuel_prices_final_1 = pd.merge(df_fuel_prices_final.loc[~(df_fuel_prices_final['TECHNOLOGY']
                                                                .str[6:9].isin(country_list))],
                                      df_int_fuel_prices,
                                      how='left',
                                      on=['FUEL', 'YEAR'])
    df_fuel_prices_final_1 = df_fuel_prices_final_1[['TECHNOLOGY',
                                                     'MODE_OF_OPERATION',
                                                     'YEAR',
                                                     'VALUE']]
    df_fuel_prices_final_1 = df_fuel_prices_final_1.pivot(index=['YEAR',
                                                                 'MODE_OF_OPERATION'],
                                                          columns='TECHNOLOGY',
                                                          values='VALUE').reset_index()
    df_fuel_prices_final_1 = df_fuel_prices_final_1.interpolate(method='linear', 
                                                                limit_direction='both')
    df_fuel_prices_final_1 = df_fuel_prices_final_1.dropna(axis='columns')
    df_fuel_prices_final_1 = pd.melt(df_fuel_prices_final_1,
                                     id_vars=['YEAR','MODE_OF_OPERATION'],
                                     value_vars=[x for x in df_fuel_prices_final_1
                                                 if x not in
                                                 ['YEAR','MODE_OF_OPERATION']],
                                     var_name='TECHNOLOGY',
                                     value_name='VALUE')
    df_fuel_prices_final_1['VALUE'] = df_fuel_prices_final_1['VALUE'].round(2)

    # Merge values for countries with country-specific data
    df_fuel_prices_final_2 = pd.merge(df_fuel_prices_final.loc[(df_fuel_prices_final['TECHNOLOGY']
                                                                .str[6:9].isin(country_list))],
                                      df_fuel_prices,
                                      how='left',
                                      on=['TECH_COUNTRY', 'YEAR'])

    df_fuel_prices_final_2 = df_fuel_prices_final_2[['TECHNOLOGY',
                                                     'MODE_OF_OPERATION',
                                                     'YEAR',
                                                     'VALUE']]
    df_fuel_prices_final_2 = df_fuel_prices_final_2.pivot(index=['YEAR',
                                                                 'MODE_OF_OPERATION'],
                                                          columns='TECHNOLOGY',
                                                          values='VALUE').reset_index()
    df_fuel_prices_final_2 = df_fuel_prices_final_2.interpolate(method='linear', 
                                                                limit_direction='both')
    df_fuel_prices_final_2 = df_fuel_prices_final_2.dropna(axis='columns')
    df_fuel_prices_final_2 = pd.melt(df_fuel_prices_final_2,
                                     id_vars=['YEAR','MODE_OF_OPERATION'],
                                     value_vars=[x for x in df_fuel_prices_final_2
                                                 if x not in
                                                 ['YEAR','MODE_OF_OPERATION']],
                                     var_name='TECHNOLOGY',
                                     value_name='VALUE')
    df_fuel_prices_final_2['VALUE'] = df_fuel_prices_final_2['VALUE'].round(2)

    # Combine dataframes with country-specific and international values
    df_fuel_prices_final = pd.concat([df_fuel_prices_final_1,
                                      df_fuel_prices_final_2])

    df_fuel_prices_final['REGION'] = region_name
    df_fuel_prices_final = df_fuel_prices_final[['REGION', 
                                                 'TECHNOLOGY', 
                                                 'MODE_OF_OPERATION',
                                                 'YEAR',
                                                 'VALUE']]
    
    return df_fuel_prices_final