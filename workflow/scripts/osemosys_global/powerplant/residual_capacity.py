"""Function to calculate residual capacity for powerplant technologies."""

import pandas as pd

from constants import (
    start_year,
    end_year,
    region_name,
    custom_nodes,
)

from data import(
    createPwrTechs, 
    custom_nodes_csv
    )


def res_capacity(df_gen_base, tech_list, df_tech_code, df_custom_res_cap, duplicate_techs):

    # ### Calculate residual capacity
    res_cap_cols = [
        "node_code",
        "tech_code",
        "total_capacity",
        "start_year",
        "retirement_year_model",
    ]

    df_res_cap = df_gen_base[res_cap_cols]
    
    df_res_cap = df_res_cap.reindex(columns=[*df_res_cap.columns.tolist(), *list(range(start_year, end_year+1))], fill_value=0)

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
    df_res_cap = createPwrTechs(df_res_cap, duplicate_techs)

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

    if custom_nodes:
        df_res_cap_custom, custom_techs = custom_nodes_csv(df_custom_res_cap, 
                                                           tech_list)
        df_res_cap = pd.concat([df_res_cap, df_res_cap_custom])

    df_res_cap.drop_duplicates(subset=['REGION','TECHNOLOGY','YEAR'],
                               keep='last',
                               inplace=True)
    df_res_cap = df_res_cap.loc[(df_res_cap['TECHNOLOGY'].str.startswith('PWR')) &
                                (~df_res_cap['TECHNOLOGY'].str.endswith('00'))]
    
    return df_res_cap, custom_techs if custom_nodes else df_res_cap