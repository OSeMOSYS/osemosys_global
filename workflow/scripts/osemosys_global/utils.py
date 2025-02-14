"""Utility Functions"""

import pandas as pd
from typing import Optional
from constants import SET_DTYPES

import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def apply_timeshift(x, timeshift):
    """Applies timeshift to organize dayparts.
    
    Arguments:
        x = Value between 0-24
        timeshift = value offset from UTC (-11 -> +12)"""

    x += timeshift
    if x > 23:
        return x - 24
    elif x < 0:
        return x + 24
    else:
        return x

def apply_dtypes(df:pd.DataFrame, name: Optional[str]) -> pd.DataFrame:
    """Sets datatypes on dataframe"""
    
    for col in df.columns:
        try:
            df[col] = df[col].astype(SET_DTYPES[col])
        except KeyError:
            if name:
                logging.info(f"Can not set dtype for {name} on {col}")
            else:
                logging.info(f"Can not set dtype on {col}")
    return df

def discount_factor(
    regions: list,
    years: list,
    discount_rate: pd.DataFrame, 
) -> pd.DataFrame:
    """DiscountFactor

    Arguments
    ---------
    regions: list
    years: list
    discount_rate: pd.DataFrame
    
    Notes
    -----
    From the formulation::

        param DiscountFactor{r in REGION, y in YEAR} :=
                (1 + DiscountRate[r]) ^ (y - min{yy in YEAR} min(yy)); 
    """

    if discount_rate.empty:
        raise ValueError(
            "Cannot calculate discount factor due to missing discount rate"
        )
 
    discount_rate["YEAR"] = [years] * len(discount_rate)
    discount_factor = discount_rate.explode("YEAR").reset_index(level="REGION")
    discount_factor["YEAR"] = discount_factor["YEAR"].astype("int64")
    discount_factor["NUM"] = discount_factor["YEAR"] - discount_factor["YEAR"].min()
    discount_factor["RATE"] = discount_factor["VALUE"] + 1
    discount_factor["VALUE"] = (
        discount_factor["RATE"].pow(discount_factor["NUM"] + adj).astype(float)
    )
    return discount_factor.reset_index()[["REGION", "YEAR", "VALUE"]].set_index(
        ["REGION", "YEAR"]
    )
  