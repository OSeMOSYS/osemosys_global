
import requests
import os
import sys
import pandas as pd
from osemosys_global.configuration import ConfigFile, ConfigPaths
from osemosys_global.utils import get_config_data
from pathlib import Path 

# LOGGING
import logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def main(input_dir: str, output_dir: str, start_year: int, end_year: int):
    """Creates capacity limits on renewable technologies.
    
    Args:
        input_data: str
            Path to resources/data folder 
        output_dir: str
            Path to results/data folder
        start_year: int
            Start year of demand 
        end_year: int
            End year of demand 
    """
    
    # set path variables
    INPUT_DATA_DIR = Path(input_dir)
    OUTPUT_DATA_DIR = Path(output_dir)
    YEARS = list(range(start_year, end_year + 1))
    REGION = "GLOBAL"

    ## Checks whether PLEXOS-World/MESSAGEix-GLOBIOM soft-link model data needs to be 
    # retrieved from the PLEXOS-World Harvard Dataverse.
    try:
        df_reslimit = pd.read_excel(os.path.join(
            INPUT_DATA_DIR, "data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx"),
            sheet_name = "Properties")

    except IOError:
        url = 'https://dataverse.harvard.edu/api/access/datafile/6040815'
        r = requests.get(url)
        with open(os.path.join(INPUT_DATA_DIR, 'data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx'), 'wb') as outfile:
            outfile.write(r.content)

        df_reslimit = pd.read_excel(os.path.join(
            INPUT_DATA_DIR, "data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx"),
            sheet_name = "Properties")

    # TECHNOLOGY MAPPING FOR PLEXOS -> OSEMOSYS GLOBAL

    dict_reslimit = {'Hydro' : 'HYD',
                     'Solar|CSP' : 'CSP',
                     'Solar|PV' : 'SPV',
                     'Wind|Onshore' : 'WON',
                     'Wind|Offshore' : 'WOF'}

    # GET CAPAPCITY LIMITS FROM PLEXOS WORLD 
    # This is capacity ADDITION limits, not total max capacity limits 

    df_reslimit_units = df_reslimit.loc[(df_reslimit['child_object'].str.contains('|'.join(dict_reslimit.keys()))) & 
                                  (df_reslimit['property'] == 'Max Units Built') &
                                  (df_reslimit['scenario'].str.contains('Base')) &
                                  (df_reslimit['child_class'] == 'Generator')].set_index('child_object')

    df_reslimit_capacity = df_reslimit.loc[(df_reslimit['child_object'].str.contains('|'.join(dict_reslimit.keys()))) & 
                                  (df_reslimit['property'] == 'Max Capacity') &
                                  (df_reslimit['child_class'] == 'Generator')].set_index('child_object')

    df_reslimit_final = pd.DataFrame(df_reslimit_capacity['value'] * df_reslimit_units['value'] / 1000).rename(columns = {'value' : 'VALUE'})
    df_reslimit_final['node'], df_reslimit_final['powerplant']  = df_reslimit_final.index.str.rsplit('|',1).str[1], df_reslimit_final.index.str.rsplit('|',1).str[0]
    df_reslimit_final['powerplant'] = df_reslimit_final['powerplant'].map(dict_reslimit)

    df_reslimit_final.loc[df_reslimit_final['node'].str.len() <= 6, 
                 'node_code'] = (df_reslimit_final['node'].
                                 str.split('-').
                                 str[1:].
                                 str.join("") +
                                 'XX')
    df_reslimit_final.loc[df_reslimit_final['node'].str.len() > 6, 
                 'node_code'] = (df_reslimit_final['node'].
                                 str.split('-').
                                 str[1:].
                                 str.join("")
                                )

    df_reslimit_final['TECHNOLOGY'] = 'PWR' + df_reslimit_final['powerplant'] + df_reslimit_final['node_code'] + '01'
    cap_addition_limit = df_reslimit_final.set_index('TECHNOLOGY').to_dict()['VALUE']

    # GET RESIDUAL CAPACITY VALUES 

    df_res_cap_raw = pd.read_csv(os.path.join(OUTPUT_DATA_DIR, 'ResidualCapacity.csv'))
    df_res_cap_raw['VALUE'] = df_res_cap_raw.loc[:,'VALUE'].round(4)
    df_res_cap = df_res_cap_raw.loc[
        df_res_cap_raw['TECHNOLOGY'].str[3:6].isin(list(dict_reslimit.values()))]
    df_res_cap = get_max_value_per_technology(df_res_cap)
    res_cap = df_res_cap.set_index('TECHNOLOGY').to_dict()['VALUE']
    
    # CALCULATE AND FORMAT DATA 

    out_data = []
    for tech, maxad in cap_addition_limit.items():
        try:
            max_capacity = res_cap[tech] + maxad
        except KeyError:
            max_capacity = maxad

        # Add 0.0002 to enusre there is no rounding mismathch between total 
        # annual max capacity and residual capacity 
        max_capacity = round(max_capacity, 4) + 0.0002

        for year in YEARS:
            out_data.append([
                REGION,
                tech,
                year,
                max_capacity
            ])

    df_max_capacity = pd.DataFrame(out_data, columns=[
        'REGION',
        'TECHNOLOGY',
        'YEAR',
        'VALUE'
    ])
    df_max_capacity.to_csv(os.path.join(
        OUTPUT_DATA_DIR, "TotalAnnualMaxCapacity.csv"), index = None)

def get_max_value_per_technology(df: pd.DataFrame) -> pd.DataFrame:
    """Gets the max value for each unique technology in a dataframe. 

    This function will search through a 'TECHNOLOGY' column to idnetify each
    unique technology. The input dataframe will be filtered based on each 
    unique technology and only keep one datapoint per technology - the 
    datapoint cooresponding to the max value in the 'VALUE' column. 

    Args: 
        df: Dataframe with at minimum a 'TECHNOLOGY' and 'VALUE' columns

    Returns: 
        df: Filtered dataframe giving max values per technology. 
    """ 

    # Get list of techs to filter over 
    techs = df['TECHNOLOGY'].unique().tolist()

    #output list to hold filtered data
    out_data = []

    # perform filtering
    for tech in techs:
    # for tech in ['PWRHYDCHNJS01']: 
        df_tech_filter = df.loc[df['TECHNOLOGY'] == tech]
        df_value_filter = df_tech_filter.loc[
            df_tech_filter['VALUE'] == df_tech_filter['VALUE'].max()].reset_index(drop=True)
        df_value_filter = df_value_filter.drop_duplicates(subset=['VALUE'])
        value_filter_data = df_value_filter.values.tolist()
        out_data.append(value_filter_data[0])
    
    # setup dataframe to return
    df_out = pd.DataFrame(out_data, columns=list(df))
    return df_out

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python demand_projection.py <resources/data> <results/data> <config.yaml>")
        logging.info("Max capacity limits failed")
    else:
        config_values = ["startYear", "endYear", "nodes_to_add"]
        config_data = get_config_data(sys.argv[3], config_values)
        
        main(
            sys.argv[1], 
            sys.argv[2], 
            int(config_data["startYear"]),
            int(config_data["endYear"]),
        )
        logging.info("Max capacity limits sucessfully set")