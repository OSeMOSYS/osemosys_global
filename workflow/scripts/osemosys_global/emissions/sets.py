"""Function to create sets."""

import pandas as pd

def set_unique_emissions(ar: pd.DataFrame) -> pd.DataFrame():
    data = ar.copy()
    data = data[["EMISSION"]].drop_duplicates()
    data.rename(columns={"EMISSION": "VALUE"}, inplace=True)
    
    return data