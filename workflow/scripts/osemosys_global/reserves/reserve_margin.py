"""Functions to set ReserveMargin."""
import pandas as pd

from data import get_years

def set_reserve_margin(reserve_margins, start_year, 
                       end_year, region_name):
    
    years = get_years(start_year, end_year)

    if reserve_margins:
        data = pd.DataFrame(years, columns=["YEAR"])
        for rm, rm_params in reserve_margins.items():
            data.loc[data["YEAR"].between(rm_params[1], rm_params[2]), "VALUE"] = (
                1 + rm_params[0] / 100
            )
    
        data = data.interpolate()
        data["REGION"] = region_name
        data = data[["REGION", "YEAR", "VALUE"]]
    else:
        data = pd.DataFrame(columns=["REGION", "YEAR", "VALUE"])
    
    return data