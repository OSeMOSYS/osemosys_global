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

def costs_var_pwr(df_cmo_data, df_fuel_prices, tech_set,
                  start_year, end_year, region_name):
    

    # ### Import data files and user input              

    ## Data for variable costs of fuels taken from World Bank Commodity Market Outlooks:
    ##    https://www.worldbank.org/en/research/commodity-markets
    ## Download the 'Charts and Data' zip file and extract the forecasts file (CMO-April-2020-forecasts.xlsx).
    ## Adjust filename below for forecasts done at a different date.

    # Read in World Bank Commodity Price Outlook - we only want rows 87 - 91
    # using 85 as the headers (years) and skipping the energy header...

    #### REPLACE FROM BELOW ONCE FIGURED OUT...df_cmo_data = pd.read_excel(

    # Read in Technologies
    years =  get_years(start_year, end_year)
    # ### Filter technologies to keep only fuel production technologies (MIN)

    df_techs = tech_set[tech_set.VALUE.str.contains('|'.join(['MIN','RNW']))]
    techs = df_techs['VALUE']

    '''
    df_techs['REGION'] = region_name
    df_techs['MODE_OF_OPERATION'] = 1
    df_techs_temp = df_techs.copy()
    df_techs_temp['MODE_OF_OPERATION'] = 2
    df_techs = pd.concat([df_techs,df_techs_temp])
    df_techs.rename(columns={'VALUE':'TECHNOLOGY'}, inplace=True)


    df_techs['YEAR'] = [range(years[0], years[-1] + 1)] * len(df_techs)

    df_techs = df_techs.explode('YEAR')

    df_techs.reset_index(drop=True, inplace=True)

    # print(df_techs)


    # ### Cleanup prices from CMO data

    df_cmo_data = df_cmo_data.drop([0], axis=0)
    df_cmo_data = df_cmo_data.drop(["Commodity"], axis=1)
    columns = df_cmo_data.columns.values
    columns[0] = "YEAR"
    df_cmo_data = df_cmo_data.loc[:, ~df_cmo_data.columns.str.contains('^Unnamed', na=False)]

    df_cmo_data = df_cmo_data.transpose()

    # use commodity key as header
    new_header = df_cmo_data.iloc[0] #grab the first row for the header
    df_cmo_data = df_cmo_data[1:] #take the data less the header row
    new_header[1] = 'MINCOA'  # This was KFCOAL_AUS  KFCRUDE_PETRO  KFNGAS_EUR
    new_header[2] = 'MINOIL'  # This was KFCRUDE_PETRO
    new_header[3] = 'MINGAS'  # This was KFNGAS_EUR
    df_cmo_data.columns = new_header #set the header row as the df header

    # drop units as we're doing $/PJ
    df_cmo_data = df_cmo_data.drop("Unit", axis = 0)

    # And convert to $/PJ
    # Values taken from kylesconverter.com:
    #  1 MT coal contains 29.31 PJ (1,000,000 tonnes coal) Original was $/mt ($/tonne I think).  We want $mill/PJ so divide by 29.31
    #  1 MMBtu = 0.000001055056 PJ Original was MMBtu and we want $mill/pJ, so divide by 1.055056.
    #  1 bbl = 0.00000612 PJ (Barrels of Oil) Original was $/bbl so divide by 0.00000612 and multiply by 1000000
    # NOT ENTIRELY SURE I'M CALCULATING THE UNITS CORRECTLY (NOT SURE WHAT ORIGINAL UNITS WERE AS IT'S NOT STATED ANYWHERE)

    # ORIGINAL UNITS:
    #       KFCOAL_AUS KFCRUDE_PETRO KFNGAS_EUR KFNGAS_US KFNGAS_JP
    #Unit       $/mt         $/bbl    $/mmbtu   $/mmbtu   $/mmbtu

    df_cmo_data['MINCOA'] = df_cmo_data['MINCOA'] / 29.31
    df_cmo_data['MINOIL'] = df_cmo_data['MINOIL'] / 0.00000612 / 1000000
    df_cmo_data['MINGAS'] = df_cmo_data['MINGAS'] / 1.055056
    df_cmo_data['KFNGAS_US'] = df_cmo_data['KFNGAS_US'] / 1.055056
    df_cmo_data['KFNGAS_JP'] = df_cmo_data['KFNGAS_JP'] / 1.055056

    # print(df_cmo_data)

    # Add in other fuels that are the same as those above:
    df_cmo_data['MINCOG'] = df_cmo_data['MINGAS']  # Cogen is powered by gas
    df_cmo_data['MINOTH'] = df_cmo_data['MINOIL']  # Other petroleum products are similar to oil
    df_cmo_data['MINPET'] = df_cmo_data['MINOIL']  # Petroleum products are similar to oil

    # Add price for URN. 40 $2020/lb -> 0.0226 m$2020/PJ (3900 GJ/kg)
    df_cmo_data['MINURN'] = 0.0226

    # And add in international prices that are 50% higher than the regular ones:
    df_cmo_data['INTCOA'] = df_cmo_data['MINCOA'] * 1.5
    df_cmo_data['INTOIL'] = df_cmo_data['MINOIL'] * 1.5
    df_cmo_data['INTGAS'] = df_cmo_data['MINGAS'] * 2
    df_cmo_data['INTCOG'] = df_cmo_data['MINCOG'] * 2
    df_cmo_data['INTOTH'] = df_cmo_data['MINOTH'] * 2
    df_cmo_data['INTPET'] = df_cmo_data['MINPET'] * 2
    df_cmo_data['INTURN'] = df_cmo_data['MINURN'] * 1.5

    df_cmo_data = df_cmo_data.reindex(years)

    df_cmo_data = df_cmo_data.apply(pd.to_numeric)

    df_cmo_data = df_cmo_data.interpolate()

    col = df_cmo_data.columns
    df_cmo_data = pd.melt(df_cmo_data.reset_index(), id_vars='index', value_vars=col)

    df_cmo_data = df_cmo_data.rename(columns={'index': 'YEAR', 'YEAR': 'TEMPTECH'})

    # print(df_cmo_data)


    # ### Map costs to technologies by region/country

    ## Need to create: REGION,TECHNOLOGY,MODE_OF_OPERATION,YEAR,VALUE

    # Setup TEMPTECH column for merge
    df_techs['TEMPTECH'] = df_techs['TECHNOLOGY'].str[0:6]

    # International costs are identified by the INT in MINCOAINT, so make these INTCOA temporarily
    df_techs.loc[df_techs.TECHNOLOGY.str[6:9]=='INT', 'TEMPTECH'] = 'INT'+df_techs['TECHNOLOGY'].str[3:6]

    df_varcost = pd.merge(df_techs, df_cmo_data, on=['YEAR', 'TEMPTECH'])

    df_varcost = df_varcost.rename(columns={'value': 'VALUE'})

    df_varcost = df_varcost.drop(['TEMPTECH'], axis=1)

    # print(df_varcost)

    ## AND THEN DROP ANY ITEMS THAT HAVE A NAN AS THESE ARE DEFAULT (0) VALUES
    ## BUT THIS IS ALREADY DONE FOR US BY THE MERGE
    #df_techs.dropna(subset = ["VALUE"], inplace=True)
    #print(df_techs)

    # ### Write out variablecost.csv

    df_varcosts_final = df_varcost[['REGION', 
                           'TECHNOLOGY', 
                           'MODE_OF_OPERATION',
                           'YEAR', 
                           'VALUE']]
    df_varcosts_final['VALUE'] = df_varcosts_final['VALUE'].round(2)
    #df_varcosts_final.to_csv(os.path.join(output_data_dir,'VariableCost.csv'), mode='w', header=True, index = None)
    '''
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

    # df_varcosts_final = apply_dtypes(df_varcosts_final, "VariableCost")
    # df_varcosts_final.to_csv(os.path.join(output_data_dir,'VariableCost.csv'), mode='w', header=True, index = None)

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