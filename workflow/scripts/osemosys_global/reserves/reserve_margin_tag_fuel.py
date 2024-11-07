"""Functions to set ReserveMarginTagFuel."""
import pandas as pd
import itertools

from data import get_years

def set_reserve_margin_fuels(fuel_set, start_year, 
                             end_year, region_name):
    
    years = get_years(start_year, end_year)

    rm_fuels = [
        x for x in fuel_set["VALUE"].astype(str).unique(
            ) if x.startswith("ELC") if x.endswith("01")
    ]
    
    data = pd.DataFrame(
        list(itertools.product([region_name], rm_fuels, years, [1])),
        columns=["REGION", "FUEL", "YEAR", "VALUE"],
    )
    
    return data