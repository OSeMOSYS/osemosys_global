"""Functions to set ReserveMarginTagTechnology."""
import pandas as pd
import itertools

from data import get_years

def set_reserve_margin_technologies(margins_technologies, 
                                    tech_set, start_year, 
                                    end_year, region_name):

    years = get_years(start_year, end_year)
    data = pd.DataFrame()

    '''Check if transmission is added to the config to contribute to
    reserves and develop outputs.'''
    for tech, val in margins_technologies.items():
        if tech == 'TRN':
            rm_techs = [
                x
                for x in tech_set["VALUE"].unique()
                if x.startswith("TRN")
                ]

        else:
            rm_techs = [
                x
                for x in tech_set["VALUE"].unique()
                if x.startswith("PWR")
                if x[3:6] == tech
            ]

        tech_data = pd.DataFrame(
            list(itertools.product([region_name], rm_techs, years)),
            columns=["REGION", "TECHNOLOGY", "YEAR"],
        )

        tech_data['VALUE'] = val / 100
        data = pd.concat([data, tech_data])

    backstop = get_backstop_rm(tech_set, region_name, start_year, end_year)

    return pd.concat([data, backstop])


def get_backstop_rm(
    tech_set: pd.DataFrame, region: str, start_year: int, end_year: int
) -> pd.DataFrame:
    """Backstop does not contribute to reserve margin"""

    bck_techs = tech_set[tech_set.VALUE.str.startswith("PWRBCK")]

    if bck_techs.empty:
        return pd.DataFrame(columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"])

    t = bck_techs.VALUE.unique().tolist()
    y = get_years(start_year, end_year)

    df = pd.DataFrame(
        index=pd.MultiIndex.from_product([t, y], names=["TECHNOLOGY", "YEAR"])
    ).reset_index()
    df["REGION"] = region
    df["VALUE"] = 0

    return df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
