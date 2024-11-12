"""Function to set EmissionsPenalty."""

import pandas as pd
import itertools

from data import get_years

def get_emission_penalty(emissions, emission_penalty, 
                         start_year, end_year, region_name):
    """Creates emission penalty dataframe.

    The emission penalty is applied at a global geographical level. All regions
    and subregions have the same penalty.

    Args:
        emission: string describing the emission type (ie. 'CO2')
        penalty: emission penalty in M$/MT

    Returns:
        df: Dataframe describing emission penalty. Dataframe headers are REGION,
            EMISSION, YEAR, VALUE
    """
    
    if emission_penalty:
        emissions_list = list(emissions["VALUE"])
        years = get_years(start_year, end_year)
    
        # GENERATE DATA
    
        # Create dataframe template to calculate SpecifiedAnnualDemand
        df = pd.DataFrame(
            list(itertools.product(emissions_list, years)), columns=["EMISSION", "YEAR"]
        )
    
        for ep_params in emission_penalty:
            df.loc[
                (df["YEAR"].between(ep_params[2], ep_params[3]))
                & (df["EMISSION"].isin([ep_params[0] + ep_params[1]])),
                "VALUE",
            ] = ep_params[4]
    
        df = df.pivot(index=["YEAR"], columns=["EMISSION"], values="VALUE").reset_index()
    
        # Drop all columns with only NaN
        df.dropna(axis=1, how="all", inplace=True)
    
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
        df = pd.DataFrame(columns=["REGION", "EMISSION", "YEAR", "VALUE"])

    return df