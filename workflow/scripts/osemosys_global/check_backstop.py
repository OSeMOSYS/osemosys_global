"""Prints warning if backstop techs are used"""

import pandas as pd
import sys
from pathlib import Path

if __name__ == "__main__":

    if len(sys.argv) != 2:
        raise ValueError("Usage: python create_checkbackstop.py <sceanrio_name>")
    else:
        scenario = sys.argv[1]

    new_capacity_csv = Path("results", scenario, "results", "NewCapacity.csv")
    try:
        new_capacity = pd.read_csv(new_capacity_csv)
    except FileNotFoundError:
        sys.exit() 

    backstop = new_capacity[new_capacity.TECHNOLOGY.str.startswith("PWRBCK")]

    if not backstop.empty:
        techs = backstop.TECHNOLOGY.unique().tolist()
        print("\n*****\n")
        print("The following Backstop Technologies are being used:\n")
        print(techs)
        print("\n*****\n")
