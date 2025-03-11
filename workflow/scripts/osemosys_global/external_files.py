import logging
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)

from pathlib import Path
import os
import requests
import sys

def download_file(file: str, url: Path | str) -> None:
    """Downloads a file if the file does not already exist"""

    if Path(file).exists():
        logging.info(f"{file} already exists")
        return

    logging.info(f"Downloading {file}")

    data = requests.get(url, path)

    with open(path, "wb") as f:
        f.write(data.content)

if __name__ == "__main__":

    if "snakemake" in globals():
        external_files = snakemake.params.files
        input_data_dir = snakemake.params.input_data_dir
    else:
        if len(sys.argv) != 3:
            msg = "Usage: python {} <save_name> <url>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            in_file = sys.argv[1]
            in_url = sys.argv[2]
            external_files = {in_file: in_url}
            
        input_data_dir = 'resources/data/default'

    for file, url in external_files.items():
        path = os.path.join(input_data_dir, file)
        download_file(file, url)