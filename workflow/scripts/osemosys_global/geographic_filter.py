# Filter osemosys_global datapackaged based on user-defined geographic scope

import pandas as pd
from typing import Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

INT_FUELS = ["COA", "COG", "GAS", "OIL", "PET", "OTH", "URN"]


def filer(
    df: pd.DataFrame,
    name: str,
    geo_scope: list[str],
    remove_nodes: list[str],
    res_targets: Optional[dict[str, list]] = None,
) -> pd.DataFrame:

    if df.empty:
        return df

    # Do not filter if only element is international fuels
    if len(geo_scope) == 1 and geo_scope[0] == "INT":
        return df

    if "TECHNOLOGY" in df.columns:
        df = df.loc[
            df["TECHNOLOGY"].str[3:6].isin(geographic_scope)
            | df["TECHNOLOGY"].str[6:9].isin(geographic_scope)
            | df["TECHNOLOGY"].str[8:11].isin(geographic_scope)
        ]

        # Filter out all international TRN techs
        df = df.loc[
            ~(
                df["TECHNOLOGY"].str.startswith("TRN")
                & (
                    ~(df["TECHNOLOGY"].str[3:6].isin(geographic_scope))
                    | ~(df["TECHNOLOGY"].str[8:11].isin(geographic_scope))
                )
            )
        ]

        if remove_nodes:
            df = df.loc[
                ~(
                    df["TECHNOLOGY"].str[3:8].isin(remove_nodes)
                    | df["TECHNOLOGY"].str[6:11].isin(remove_nodes)
                    | df["TECHNOLOGY"].str[8:13].isin(remove_nodes)
                )
            ]

    if "STORAGE" in df.columns:
        df = df.loc[
            df["STORAGE"].str[3:6].isin(geographic_scope)
            | df["STORAGE"].str[6:9].isin(geographic_scope)
            | df["STORAGE"].str[8:11].isin(geographic_scope)
        ]

        if remove_nodes:
            df = df.loc[
                ~(
                    df["STORAGE"].str[3:8].isin(remove_nodes)
                    | df["STORAGE"].str[6:11].isin(remove_nodes)
                    | df["STORAGE"].str[8:13].isin(remove_nodes)
                )
            ]

    if "FUEL" in df.columns:
        if res_targets is None:
            df = df.loc[
                df["FUEL"].str[3:6].isin(geographic_scope)
                | df["FUEL"].str[6:9].isin(geographic_scope)
                | df["FUEL"].isin(INT_FUELS)
            ]

        else:
            df = df.loc[
                df["FUEL"].str[3:6].isin(geographic_scope)
                | df["FUEL"].str[6:9].isin(geographic_scope)
                | df["FUEL"].isin(res_targets)
                | df["FUEL"].isin(INT_FUELS)
            ]

        if remove_nodes:
            df = df.loc[
                ~(
                    df["FUEL"].str[3:8].isin(remove_nodes)
                    | df["FUEL"].str[6:11].isin(remove_nodes)
                )
            ]

    if name == "FUEL":
        if res_targets is None:
            df = df.loc[
                df["VALUE"].str[3:6].isin(geographic_scope)
                | df["VALUE"].str[6:9].isin(geographic_scope)
                | df["VALUE"].isin(INT_FUELS)
            ]

        else:
            df = df.loc[
                df["VALUE"].str[3:6].isin(geographic_scope)
                | df["VALUE"].str[6:9].isin(geographic_scope)
                | df["VALUE"].isin(res_targets)
                | df["VALUE"].isin(INT_FUELS)
            ]

        if remove_nodes:
            df = df.loc[
                ~(
                    df["VALUE"].str[3:8].isin(remove_nodes)
                    | df["VALUE"].str[6:11].isin(remove_nodes)
                )
            ]

    if name == "TECHNOLOGY":
        df = df.loc[
            df["VALUE"].str[3:6].isin(geographic_scope)
            | df["VALUE"].str[6:9].isin(geographic_scope)
            | df["VALUE"].str[8:11].isin(geographic_scope)
        ]
        df = df.loc[
            ~(
                df["VALUE"].str.startswith("TRN")
                & (
                    ~(df["VALUE"].str[3:6].isin(geographic_scope))
                    | ~(df["VALUE"].str[8:11].isin(geographic_scope))
                )
            )
        ]

        if remove_nodes:
            df = df.loc[
                ~(
                    df["VALUE"].str[3:8].isin(remove_nodes)
                    | df["VALUE"].str[6:11].isin(remove_nodes)
                    | df["VALUE"].str[8:13].isin(remove_nodes)
                )
            ]

    if name == "STORAGE":
        df = df.loc[df["VALUE"].str[3:6].isin(geographic_scope)]

        if remove_nodes:
            df = df.loc[~df["VALUE"].str[3:8].isin(remove_nodes)]

    return df


if __name__ == "__main__":

    if "snakemake" in globals():
        geographic_scope = snakemake.params.geographic_scope
        res_targets = snakemake.params.res_targets
        nodes_to_remove = snakemake.params.nodes_to_remove
        in_dir = snakemake.params.in_dir
        out_dir = snakemake.params.out_dir
    else:
        geographic_scope = ["IND"]
        res_targets = {"T01": ["", [], "PCT", 2048, 2050, 95]}
        nodes_to_remove = []
        in_dir = "results/data"
        out_dir = "results/data/Wrong/data"

    geographic_scope.append("INT")  # for international fuels added by default

    if not Path(out_dir).exists():
        Path(out_dir).mkdir(parents=True)

    for each_csv in Path(in_dir).glob("*.csv"):
        df = pd.read_csv(Path(each_csv))
        df = filer(df, each_csv.stem, geographic_scope, nodes_to_remove, res_targets)
        df.to_csv(Path(out_dir, each_csv.name), index=False)

    logging.info("Geographic Filter Applied")
