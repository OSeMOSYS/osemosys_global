"""Function to set EmissionActivityRatio."""

import pandas as pd

from data import get_co2_emission_factors

def get_ear(emission, emission_factors, ccs_efficiency,
            iar_base, oar_base, tech_to_fuel):
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

    # GET EMISSION FACTORS

    co2_factors = get_co2_emission_factors(emission_factors)

    # GET INFO FROM OUTPUT ACTIVITY RATIO
    df = oar_base.drop(["FUEL", "VALUE"], axis=1)
    df = df[
        (df["TECHNOLOGY"].str.startswith("PWR"))
        & ~(df["TECHNOLOGY"].str.startswith("PWRTRN"))
    ]

    # ADD MAPPING OF TECHNOLOGY TO EMISSION ACTIVITY RATIO

    df["TECH_CODE"] = df["TECHNOLOGY"].str[3:6]
    df["COUNTRY"] = df["TECHNOLOGY"].str[6:9]
    df["FUEL_NAME"] = df["TECH_CODE"].map(tech_to_fuel)
    df["VALUE"] = df["FUEL_NAME"].map(co2_factors)
    
    # Drop all non-CO2 emitting techs
    df = df[~df.FUEL_NAME.isna()]
    
    # Multiply by InputActivityRatio
    df_iar = iar_base.copy()
    df_iar.rename(columns={"VALUE": "IAR"}, inplace=True)
    df = pd.merge(
        df, df_iar, how="left", on=["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"]
    )

    df["VALUE"] = df["VALUE"] * df["IAR"]
    df.drop_duplicates(inplace=True)
    df["VALUE"] = df["VALUE"].round(4)

    # ADD IN EMISSION COLUMN

    df["EMISSION"] = emission + df["COUNTRY"]
    
    # Add CCS EAR as a function of COA EAR and efficiency constant
    df.loc[df["TECH_CODE"].str.startswith("CCS"), "VALUE"] = (
        df.loc[df["TECH_CODE"].str.startswith("COA"), "VALUE"].mean() * (
            (100 - ccs_efficiency) / 100)
    )
    df["VALUE"] = df["VALUE"].round(4)
    # Final EmissionActivityRatio dataframe
    df = df[["REGION", "TECHNOLOGY", "EMISSION", "MODE_OF_OPERATION", "YEAR", "VALUE"]]
    
    return df