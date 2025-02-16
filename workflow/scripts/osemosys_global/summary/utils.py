import pandas as pd
from pathlib import Path


def get_discount_factor(
    years: list,
    discount_rate: pd.DataFrame,
) -> pd.DataFrame:
    """DiscountFactor

    Arguments
    ---------
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
        discount_factor["RATE"].pow(discount_factor["NUM"]).astype(float)
    )
    return discount_factor.reset_index()[["REGION", "YEAR", "VALUE"]].set_index(
        ["REGION", "YEAR"]
    )


def _read_result_data(result_dir: str) -> dict[str, pd.DataFrame]:
    """Reads in result CSVs

    Issue with otoole config
    """
    data = {}
    files = [Path(x) for x in Path(result_dir).iterdir()]
    for f in files:
        df = pd.read_csv(f)
        df = df.set_index(df.columns[:-1].tolist())
        data[f.stem] = df
    return data
