"""Sets backstop parameters"""

import pandas as pd


def get_backstop_data(
    tech_set: pd.DataFrame, year_set: pd.DataFrame, region: str
) -> tuple[pd.DataFrame]:
    """Gets the following backstop data

    - Technologies
    - OAR
    - Capital Cost
    - Fixed Cost
    - CapacityToActivity
    """

    df = tech_set.copy()

    # technologies
    techs = df[df.VALUE.str.startswith("PWR")].copy()  # pwrtrn not added yet
    techs["VALUE"] = "BCK" + df.VALUE.str[-7:-2]
    techs = techs.drop_duplicates()
    bck_techs = techs.VALUE.to_list()

    # activity ratios
    years = year_set.VALUE.to_list()
    oar = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [bck_techs, years], names=["TECHNOLOGY", "YEAR"]
        )
    ).reset_index()
    oar["REGION"] = region
    oar["FUEL"] = oar.TECHNOLOGY.str.replace("BCK", "ELC")
    oar["FUEL"] = oar.FUEL + "02"
    oar["MODE_OF_OPERATION"] = 1
    oar["VALUE"] = 1
    oar = oar[["REGION", "TECHNOLOGY", "FUEL", "MODE_OF_OPERATION", "YEAR", "VALUE"]]

    # capital cost
    capex = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [bck_techs, years], names=["TECHNOLOGY", "YEAR"]
        )
    ).reset_index()
    capex["REGION"] = region
    capex["VALUE"] = 999999
    capex = capex[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]

    # fixed costs
    opex = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [bck_techs, years], names=["TECHNOLOGY", "YEAR"]
        )
    ).reset_index()
    opex["REGION"] = region
    opex["VALUE"] = 999999
    opex = opex[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]

    # capacity to activity
    capact = pd.DataFrame(index=pd.Index(bck_techs, name="TECHNOLOGY")).reset_index()
    capact["REGION"] = region
    capact["VALUE"] = 31.536
    capact = capact[["REGION", "TECHNOLOGY", "VALUE"]]

    return techs, oar, capex, opex, capact
