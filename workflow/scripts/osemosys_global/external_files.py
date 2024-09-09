# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:14:43 2024

@author: maart
"""

import logging

logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)
from pathlib import Path
from configuration import ConfigPaths
import os
import requests
import sys


def download_file(file: str, url: Path | str) -> None:
    """Downloads a file if the file does not already exist"""

    if Path(path).exists():
        logging.info(f"{file} already exists")
        return

    logging.info(f"Downloading {file}")

    data = requests.get(url, path)

    with open(path, "wb") as f:
        f.write(data.content)


# CONFIGURATION PARAMETERS
config_paths = ConfigPaths()
input_data_dir = config_paths.input_data_dir

if __name__ == "__main__":

    if "snakemake" in globals():
        external_files = snakemake.params.files
    else:
        if len(sys.argv) != 3:
            msg = "Usage: python {} <save_name> <url>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            in_file = sys.argv[1]
            in_url = sys.argv[2]
            external_files = {in_file: in_url}

    for file, url in external_files.items():
        path = os.path.join(input_data_dir, file)
        download_file(file, url)
