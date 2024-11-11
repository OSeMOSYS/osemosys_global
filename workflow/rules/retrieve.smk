
EXTERNAL_FILES = [
    "PLEXOS_World_2015_Gold_V1.1.xlsx",
    "All_Demand_UTC_2015.csv",
    "PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx",
    "ember_yearly_electricity_data.csv",
    "GTD_existing.csv",
    "GTD_planned.csv",
    "CMO-October-2024-Forecasts.xlsx"
]

def get_external_links() -> dict[str,str]:
    """Gets links that can be downloaded via requests"""

    return {
        "PLEXOS_World_2015_Gold_V1.1.xlsx" : 
        "https://dataverse.harvard.edu/api/access/datafile/4008393?format=original&gbrecs=true",

        "All_Demand_UTC_2015.csv" :
        "https://dataverse.harvard.edu/api/access/datafile/3985039?format=original&gbrecs=true",

        "PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx" :
        "https://dataverse.harvard.edu/api/access/datafile/6040815",

        "ember_yearly_electricity_data.csv" :
        "https://storage.googleapis.com/emb-prod-bkt-publicdata/public-downloads/yearly_full_release_long_format.csv",
        
        "GTD_existing.csv" :
        "https://zenodo.org/records/10870602/files/GTD-v1.1_regional_existing.csv?download=1",
        
        "GTD_planned.csv" :
        "https://zenodo.org/records/10870602/files/GTD-v1.1_regional_planned.csv?download=1",

        "CMO-October-2024-Forecasts.xlsx" : 
        "https://thedocs.worldbank.org/en/doc/1c533d236291c53757b147f8acb8a82f-0050012024/related/CMO-October-2024-Forecasts.xlsx"
    }

rule download_external_files:
    message:
        "Downloading external files..."
    params:
        files = get_external_links()
    log:
        log = "results/logs/external_files.log"
    output:
        csv_files = expand("resources/data/{output_file}", output_file=EXTERNAL_FILES),
    script:
        "../scripts/osemosys_global/external_files.py"