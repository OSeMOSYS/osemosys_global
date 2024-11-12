"""Functions to extract and format relevent data."""

import pandas as pd

def get_years(start: int, end: int) -> range:
    return range(start, end + 1)

def _format_ember_emission_data(ember: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure

    No unit conversion needed, as Ember emissions in mtCO2
    """

    df = ember.copy()

    df = df[
        (df.Category == "Power sector emissions") & (df.Subcategory == "Total")
    ].copy()
    df["EMISSION"] = df.COUNTRY
    df["REGION"] = "GLOBAL"
    df = df[["REGION", "EMISSION", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "EMISSION", "YEAR"]).sum()

def get_co2_emission_factors(emission_factors):
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
    # https://www.epa.gov/sites/default/files/2018-03/douments/emission-factors_mar_2018_0.pdf

    # Read in emission factors
    df_raw = emission_factors.copy()
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