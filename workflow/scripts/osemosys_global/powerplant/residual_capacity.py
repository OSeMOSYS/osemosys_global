"""Function to calculate residual capacity for powerplant technologies."""

import pandas as pd
import itertools

from data import(
    create_pwr_techs, 
    get_years
    )

def res_capacity(df_gen_base, duplicate_techs, start_year, end_year, region_name):

    # ### Calculate residual capacity
    res_cap_cols = [
        "node_code",
        "tech_code",
        "total_capacity",
        "start_year",
        "retirement_year_model",
    ]

    df_res_cap = df_gen_base[res_cap_cols]
    
    df_res_cap = df_res_cap.reindex(columns=
                                    [*df_res_cap.columns.tolist(), 
                                     *list(range(start_year, end_year+1))], 
                                    fill_value=0)

    df_res_cap = pd.melt(
        df_res_cap,
        id_vars=res_cap_cols,
        value_vars=[x for x in df_res_cap.columns if x not in res_cap_cols],
        var_name="model_year",
        value_name="value",
    )
    
    df_res_cap["model_year"] = df_res_cap["model_year"].astype(int)
    df_res_cap["value"] = df_res_cap["value"].astype(float)
    
    df_res_cap.loc[
        (df_res_cap["model_year"] >= df_res_cap["start_year"])
        & (df_res_cap["model_year"] <= df_res_cap["retirement_year_model"]),
        "value",
    ] = df_res_cap["total_capacity"]

    df_res_cap = df_res_cap.groupby(
        ["node_code", "tech_code", "model_year"], as_index=False
    )["value"].sum()

    # Add column with naming convention
    df_res_cap = create_pwr_techs(df_res_cap, duplicate_techs)

    # Convert total capacity from MW to GW
    df_res_cap['value'] = df_res_cap['value'].div(1000)

    # Rename 'model_year' to 'year' and 'total_capacity' to 'value'
    df_res_cap.rename({'model_year': 'YEAR',
                       'value': 'VALUE'},
                      inplace=True,
                      axis=1)
    # Drop 'tech_code' and 'node_code'
    df_res_cap.drop(['tech_code', 'node_code'], 
                    inplace=True, 
                    axis=1)

    # Add 'REGION' column and fill 'GLOBAL' throughout
    df_res_cap['REGION'] = region_name

    # Reorder columns
    df_res_cap = df_res_cap[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]

    return df_res_cap

def add_custom_res_cap(df_res_cap, df_custom, tech_list, custom_nodes,
                       start_year, end_year, region_name):

    df_custom_res_cap = pd.DataFrame(list(itertools.product(custom_nodes,
                                                   tech_list,
                                                   get_years(start_year, end_year))
                                  ),
                             columns = ['CUSTOM_NODE',
                                        'FUEL_TYPE',
                                        'YEAR']
                             )
    df_custom_res_cap['REGION'] = region_name
    df_custom = df_custom.groupby(['CUSTOM_NODE',
                                   'FUEL_TYPE',
                                   'START_YEAR',
                                   'END_YEAR'],
                                  as_index=False)['CAPACITY'].sum()
    df_custom_res_cap = pd.merge(df_custom_res_cap,
                        df_custom,
                        how='left',
                        on=['CUSTOM_NODE',
                            'FUEL_TYPE'])
    df_custom_res_cap['TECHNOLOGY'] = ('PWR' +
                              df_custom_res_cap['FUEL_TYPE'] + 
                              df_custom_res_cap['CUSTOM_NODE'] +
                              '01')
    technologies = df_custom_res_cap['TECHNOLOGY'].unique()
    df_custom_res_cap.dropna(inplace=True)
    df_custom_res_cap.drop_duplicates(inplace=True)
    df_custom_res_cap = df_custom_res_cap.loc[df_custom_res_cap['YEAR'] >=
                                              df_custom_res_cap['START_YEAR']]
    df_custom_res_cap = df_custom_res_cap.loc[df_custom_res_cap['YEAR'] <=
                                              df_custom_res_cap['END_YEAR']]
    df_custom_res_cap['VALUE'] = df_custom_res_cap['CAPACITY'].div(1000)
    df_custom_res_cap['REGION'] = region_name
    df_custom_res_cap = df_custom_res_cap[['REGION','TECHNOLOGY','YEAR','VALUE']]
    df_custom_res_cap = df_custom_res_cap.groupby(['REGION',
                                 'TECHNOLOGY',
                                 'YEAR'],
                                 as_index=False)['VALUE'].sum()
    
    df_res_cap = pd.concat([df_res_cap, df_custom_res_cap])
    
    df_res_cap.drop_duplicates(subset=['REGION','TECHNOLOGY','YEAR'],
                               keep='last',
                               inplace=True)
    df_res_cap = df_res_cap.loc[(df_res_cap['TECHNOLOGY'].str.startswith('PWR')) &
                                (~df_res_cap['TECHNOLOGY'].str.endswith('00'))]
    
    return df_custom_res_cap, technologies