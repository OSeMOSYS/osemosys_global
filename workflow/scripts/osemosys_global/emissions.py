"""Generates emission set, emission activity ratio, and emission penalty"""

import logging
import pandas as pd
from pathlib import Path
from configuration import ConfigFile, ConfigPaths
import itertools


# Logging formatting
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)

# Constant for tech to fuel emission type mapping
_TECH_TO_FUEL = {
    #'BIO':'Bagasse',
    "WAS": "Municipal Solid Waste",
    "COA": "Lignite Coal",
    "COG": "Lignite Coal",
    "OCG": "Natural Gas",
    "CCG": "Natural Gas",
    "GAS": "Natural Gas",
    "PET": "Crude Oil",
    "OIL": "Crude Oil",
    "OTH": "Natural Gas",
    "CCS": "Lignite Coal",
}

# Emission name
_EMISSION = "CO2"


def main():
    """Assigns CO2 equivalent values to each technology over all its modes
    of operation. For technologies that do not have emissions, it assigns a
    value of zero. A global emission penalty value is applied and the
    emission type (CO2) is written to the emission set.
    """

    # CONFIGURATION PARAMETERS

    config_paths = ConfigPaths()
    config = ConfigFile("config")

    output_data_dir = config_paths.output_data_dir
    emission_penalty = config.get("emission_penalty")  # M$/MT
    emission_limit = config.get("emission_limit")  # MT of CO2-eq.
    storage_techs = config.get("storage_parameters").keys()
    storage_techs = ['PWR' + stor for stor in storage_techs]
    
    # ASSIGN EMISSION ACTIVITY RATIOS

    df_ear = get_ear(_EMISSION, storage_techs)
    df_ear.to_csv(Path(output_data_dir, "EmissionActivityRatio.csv"), index=False)
    logging.info("Successfully generated emission activity ratio")

    # ASSIGN EMISSION

    df_emission = df_ear[["EMISSION"]].drop_duplicates()
    df_emission.rename(columns={"EMISSION": "VALUE"}, inplace=True)
    df_emission.to_csv(Path(output_data_dir, "EMISSION.csv"), index=False)
    logging.info("Successfully generated emission set")

    # Create of list of EMISSIONS
    emissions = list(df_emission["VALUE"])

    # ASSIGN EMISSION PENALTY

    if emission_penalty:
        df_emission_penalty = get_emission_penalty(emissions, emission_penalty)
    else:
        df_emission_penalty = pd.DataFrame(
            columns=["REGION", "EMISSION", "YEAR", "VALUE"]
        )
    df_emission_penalty.to_csv(
        Path(output_data_dir, "EmissionsPenalty.csv"), index=False
    )
    logging.info("Successfully generated emission penalty")

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


def get_co2_emission_factors():
    """Gets co2 emission factors for diferent fuels.

    Reads in a file containing co2, ch4, n2o emission factors and global
    warming potentials for various fuel types. This function performs unit
    conversions to convert everyting to MegaTonnes per PetaJoule and collapses
    the different emission types to a single co2 equivalent value for each
    fuel.

    Returns:
        Dictionary holding fuel type as key and co2 factor as value in units
        of MT/PJ

    Example:
        co2_factors = get_co2_emission_factors()
        co2_factors['Natural Gas']
        -> 0.0503
    """

    # Source for all emission factors comes from:
    # https://www.epa.gov/sites/default/files/2018-03/documents/emission-factors_mar_2018_0.pdf

    # Configuration parameters
    config_paths = ConfigPaths()

    # Read in emission factors
    input_data_dir = config_paths.input_data_dir
    df_raw = pd.read_csv(Path(input_data_dir, "emission_factors.csv"))
    df_raw = df_raw.drop([0]).reset_index(drop=True)  # drop units row

    # Convert co2 factors from kg/mmbtu to MT/PJ
    # kg/mmbtu * 1mmbtu/1.05GJ * 1000000GJ / PJ * 1T/1000kg * 1MT/1000000T
    # Multiply by global warming potential to get co2_eq
    co2 = (
        df_raw["co2_factor"].astype(float)
        * (1 / 1055)
        * df_raw["co2_gwp"].astype(float)
    )

    # Convert ch4 and n2o factors from g/mmbtu to MT/PJ
    # kg/mmbtu * 1mmbtu/1.05GJ * 1000000GJ / PJ * 1T/1000000g * 1MT/1000000T
    # Multiply by global warming potential to get co2_eq
    ch4 = (
        df_raw["ch4_factor"].astype(float)
        * (1 / 1055000)
        * df_raw["ch4_gwp"].astype(float)
    )
    n2o = (
        df_raw["n2o_factor"].astype(float)
        * (1 / 1055000)
        * df_raw["n2o_gwp"].astype(float)
    )

    # Find total CO2 equivalent
    data = {"co2": co2, "ch4": ch4, "n2o": n2o}
    df = pd.DataFrame(data).set_axis(df_raw["FUEL TYPE"])
    df["co2_eq"] = round(df.sum(axis=1), 4).astype(float)

    return df.set_index(df.index).to_dict()["co2_eq"]


def get_ear(emission, storage_techs):
    """Creates emission activity ratio dataframe.

    This function reads in an existing input activity ratio parameter file
    and removes the fuel and year columns. This leaves a dataframe with info
    on when all technologies are allowed to operate over the model horizon.
    A column is added in to hold the emission type and emission activity ratio
    based.

    Args:
        emission: string describing the emission type (ie. 'CO2')

    Returns:
        df: Dataframe describing emission activity ratio. Dataframe headers
            are REGION, TECHNOLOGY, EMISSION, MODE_OF_OPERATION, YEAR, VALUE
    """

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    output_data_dir = config_paths.output_data_dir

    # GET EMISSION FACTORS

    co2_factors = get_co2_emission_factors()

    # GET INFO FROM INPUT ACTIVITY RATIO

    df_oar = pd.read_csv(Path(output_data_dir, "OutputActivityRatio.csv"))
    df = df_oar.drop(["FUEL", "VALUE"], axis=1)
    df = df[
        (df["TECHNOLOGY"].str.startswith("PWR"))
        & ~(df["TECHNOLOGY"].str.startswith("PWRTRN"))
    ]

    # ADD MAPPING OF TECHNOLOGY TO EMISSION ACTIVITY RATIO

    df["TECH_CODE"] = df["TECHNOLOGY"].str[3:6]
    df["COUNTRY"] = df["TECHNOLOGY"].str[6:9]
    df["FUEL_NAME"] = df["TECH_CODE"].map(_TECH_TO_FUEL)
    df["VALUE"] = df["FUEL_NAME"].map(co2_factors)
    """
    ccs_co2_factor = df.loc[df['TECH_CODE'].str.startswith('COA'),
                            'VALUE'].mean()
    ccs_co2_factor = round(ccs_co2_factor*(-3), 4)
    df.loc[df['TECH_CODE'].str.startswith('CCS'),
           'VALUE'] = ccs_co2_factor
    """
    
    # Multiply by InputActivityRatio
    df_iar = pd.read_csv(Path(output_data_dir, "InputActivityRatio.csv"))
    df_iar.rename(columns={"VALUE": "IAR"}, inplace=True)
    df = pd.merge(
        df, df_iar, how="left", on=["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"]
    )
    df["VALUE"] = df["VALUE"].fillna(0)
    df["VALUE"] = df["VALUE"] * df["IAR"]
    df.drop_duplicates(inplace=True)
    df["VALUE"] = df["VALUE"].round(4)

    # ADD IN EMISSION COLUMN

    df["EMISSION"] = emission + df["COUNTRY"]

    df.loc[df["TECH_CODE"].str.startswith("CCS"), "VALUE"] = (
        df.loc[df["TECH_CODE"].str.startswith("COA"), "VALUE"].mean() * 0.1
    )
    df["VALUE"] = df["VALUE"].round(4)
    # Final EmissionActivityRatio dataframe
    df = df[["REGION", "TECHNOLOGY", "EMISSION", "MODE_OF_OPERATION", "YEAR", "VALUE"]]
    
    # Filter out storage techs
    df = df[~df["TECHNOLOGY"].str.startswith(tuple(storage_techs))]
    
    return df


def get_emission_penalty(emissions, emission_penalty):
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
    df["REGION"] = region
    df = df[["REGION", "EMISSION", "YEAR", "VALUE"]]

    return df

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
            (df["YEAR"] == el_params[2])
            & (df["EMISSION"].isin([el_params[0] + el_params[1]])),
            "VALUE",
        ] = el_params[3]

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
