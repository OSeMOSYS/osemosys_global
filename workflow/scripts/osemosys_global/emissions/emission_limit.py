"""Function to set AnnualEmissionLimit."""

import pandas as pd
import itertools

from data import get_years

def add_emission_limits(emissions, emission_limit, start_year, 
                        end_year, region_name):

    if emission_limit:
        emissions_list = list(emissions["VALUE"])
        years = get_years(start_year, end_year)

        # GENERATE DATA
    
        # Create dataframe template
        df = pd.DataFrame(
            list(itertools.product(emissions_list, years)), columns=["EMISSION", "YEAR"]
        )
    
        for el_params in emission_limit:
            df.loc[
                (df["YEAR"] == el_params[3])
                & (df["EMISSION"].isin([el_params[0] + el_params[1]])),
                "VALUE",
            ] = el_params[4]
    
        df = df.pivot(index=["YEAR"], columns=["EMISSION"], values="VALUE").reset_index()

        # Drop all columns with only NaN
        df.dropna(axis=1, how="all", inplace=True)
        df = df.interpolate()
        # Melt to get 'COUNTRY' column and remove all rows with NaN
        df = pd.melt(
            df,
            id_vars=["YEAR"],
            value_vars=[x for x in df.columns if x not in ["YEAR"]],
            var_name="EMISSION",
            value_name="VALUE",
        )
    
        df.dropna(axis=0, inplace=True)
        df["REGION"] = region_name
        
        df = df[["REGION", "EMISSION", "YEAR", "VALUE"]]

    else:
        
        df = pd.DataFrame(
            columns=["REGION", "EMISSION", "YEAR", "VALUE"]
        )
    return df

"""

el_years = {}
years = list(range(start_year, 
                   end_year+1))
if not emission_limit is None:
    for el_params in emission_limit:
        el_years[el_params[0]] = el_params[1]

    if len(el_years) > 1:
        el_years_max =  max(el_years)
    else:
        el_years_max = int(end_year) 

    df = pd.DataFrame(list(range(min(el_years),
                                el_years_max + 1)),
                    columns=['YEAR'])
    df.sort_values(by=['YEAR'],
                inplace=True)
    df['VALUE'] = df['YEAR'].map(el_years)
    df['VALUE'].interpolate(inplace=True)
    df['VALUE'] = df['VALUE'].round(0)
    df = df[df['YEAR'].isin(years)]
    
    df['EMISSION'] = emission
    df['REGION'] = region

    df = df[['REGION',
            'EMISSION',
            'YEAR',
            'VALUE']]
else:
    df = pd.DataFrame(columns=['REGION',
                               'EMISSION',
                               'YEAR',
                               'VALUE'])
"""