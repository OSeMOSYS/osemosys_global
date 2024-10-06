"""Creates empty CSVs that are required for otoole processing"""

import yaml
import pandas as pd
from typing import Any
from pathlib import Path


def get_otoole_params(otoole_config: str) -> dict[str, dict[str, Any]]:
    """Gets parameter/result files to be created"""

    with open(otoole_config) as f:
        otoole = yaml.safe_load(f)

    return {x: otoole[x] for x in otoole if otoole[x]["type"] == "param"}


def get_empty_df(data: dict[str, Any]) -> pd.DataFrame:
    """Gets empty parameter dataframe with correct indices"""

    indices = data["indices"] + "VALUE"
    return pd.DataFrame(columns=indices)


if __name__ == "__main__":
    if "snakemake" in globals():
        otoole_yaml = snakemake.input.otoole_config
        out_dir = snakemake.params.out_dir
    else:
        otoole_yaml = "../../../resources/otoole/config.yaml"
        out_dir = "results/data"

    parameters = get_otoole_params(otoole_yaml)

    for param, data in parameters.items():
        p = Path(out_dir, f"{param}.csv")
        if not p.exists():
            df = get_empty_df(data)
            df.to_csv(str(p), index=False)
