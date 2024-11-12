"""Generates emission set, emission activity ratio, and emission penalty"""

import logging
import pandas as pd
from pathlib import Path
from configuration import ConfigFile, ConfigPaths
import itertools


# Logging formatting
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)


def main():
    """Assigns CO2 equivalent values to each technology over all its modes
    of operation. For technologies that do not have emissions, it assigns a
    value of zero. A global emission penalty value is applied and the
    emission type (CO2) is written to the emission set.
    """

    # CONFIGURATION PARAMETERS


    # ASSIGN EMISSION

    df_emission = df_ear[["EMISSION"]].drop_duplicates()
    df_emission.rename(columns={"EMISSION": "VALUE"}, inplace=True)
    df_emission.to_csv(Path(output_data_dir, "EMISSION.csv"), index=False)
    logging.info("Successfully generated emission set")

    # Create of list of EMISSIONS
    emissions = list(df_emission["VALUE"])

    # ADD EMISSION LIMITS
    if emission_limit:
        df_emission_limits = add_emission_limits(emissions, emission_limit)
    else:
        df_emission_limits = pd.DataFrame(
            columns=["REGION", "EMISSION", "YEAR", "VALUE"]
        )
    df_emission_limits.to_csv(
        Path(output_data_dir, "AnnualEmissionLimit.csv"), index=False
    )
    logging.info("Successfully generated annual emissions limit")

def add_emission_limits(emissions, emission_limit):

    # CONFIGURATION PARAMETERS
    config = ConfigFile("config")
    start_year = config.get("startYear")
    end_year = config.get("endYear")
    years = list(range(start_year, end_year + 1))
    region = config.region_name

    # GENERATE DATA

    # Create dataframe template to calculate SpecifiedAnnualDemand
    df = pd.DataFrame(
        list(itertools.product(emissions, years)), columns=["EMISSION", "YEAR"]
    )

    for el_params in emission_limit:
        df.loc[
            (df["YEAR"] == el_params[4])
            & (df["EMISSION"].isin([el_params[0] + el_params[1]])),
            "VALUE",
        ] = el_params[5]

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
    df["REGION"] = region
    df = df[["REGION", "EMISSION", "YEAR", "VALUE"]]
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
    return df


if __name__ == "__main__":
    main()
