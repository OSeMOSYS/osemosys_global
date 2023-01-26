
import requests
import os
import yaml
import pandas as pd
from OPG_configuration import ConfigFile, ConfigPaths
import itertools

# LOGGING
import logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def main():
    """Creates capacity limits on renewable technologies.
    """
    
    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile('config')  

    input_dir = config_paths.input_dir
    output_data_dir = config_paths.output_data_dir
    custom_nodes_dir = config_paths.custom_nodes_dir
    region = config.region_name
    years = config.get_years()
    custom_nodes = config.get('nodes_to_add')
    max_build = config.get('powerplant_build_rates')

    ## Checks whether PLEXOS-World/MESSAGEix-GLOBIOM soft-link model data needs to be 
    # retrieved from the PLEXOS-World Harvard Dataverse.
    try:
        df_reslimit = pd.read_excel(os.path.join(
            input_dir, "data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx"),
            sheet_name = "Properties")

    except IOError:
        url = 'https://dataverse.harvard.edu/api/access/datafile/6040815'
        r = requests.get(url)
        with open(os.path.join(input_dir, 'data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx'), 'wb') as outfile:
            outfile.write(r.content)

        df_reslimit = pd.read_excel(os.path.join(
            input_dir, "data/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx"),
            sheet_name = "Properties")

    # TECHNOLOGY MAPPING FOR PLEXOS -> OSEMOSYS GLOBAL

    dict_reslimit = {'Hydro' : 'HYD',
                     'Solar|CSP' : 'CSP',
                     'Solar|PV' : 'SPV',
                     'Wind|Onshore' : 'WON',
                     'Wind|Offshore' : 'WOF'}

    # GET CAPACITY LIMITS FROM PLEXOS WORLD 
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

    if custom_nodes:
        df_re_potentials_custom = pd.read_csv(os.path.join(custom_nodes_dir,
                                              'RE_potentials.csv')
                                              )
        df_re_potentials_custom['TECHNOLOGY'] = ('PWR' + 
                                                 df_re_potentials_custom['FUEL_TYPE'] + 
                                                 df_re_potentials_custom['CUSTOM_NODE'] + 
                                                 '01')
        re_potentials_dict = dict(zip(list(df_re_potentials_custom['TECHNOLOGY']),
                                      list(df_re_potentials_custom['CAPACITY'])
                                      )
                                  )
        cap_addition_limit.update(re_potentials_dict)

    # GET RESIDUAL CAPACITY VALUES 

    df_res_cap_raw = pd.read_csv(os.path.join(output_data_dir, 'ResidualCapacity.csv'))
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

        for year in years:
            out_data.append([
                region,
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
        output_data_dir, "TotalAnnualMaxCapacity.csv"), index = None)
    
    apply_build_rates(region, years, output_data_dir, max_build)


def get_max_value_per_technology(df):
    '''Gets the max value for each unique technology in a dataframe. 

    This function will search through a 'TECHNOLOGY' column to idnetify each
    unique technology. The input dataframe will be filtered based on each 
    unique technology and only keep one datapoint per technology - the 
    datapoint cooresponding to the max value in the 'VALUE' column. 

    Args: 
        df: Dataframe with at minimum a 'TECHNOLOGY' and 'VALUE' columns

    Returns: 
        df: Filtered dataframe giving max values per technology. 
    ''' 

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

def apply_build_rates(region, years, output_data_dir, max_build):
    """To add description

    Args:
        region:
        years:
        output_data_dir:
        max_build:

    Returns
        None
    """
    max_build_list = []
    
    if not max_build is None:
        
        max_build_df = pd.DataFrame(columns=['TYPE', 
                                             'METHOD', 
                                             'MAX_BUILD',
                                             'YEAR'])
        for tech, tech_params in max_build.items():
            '''
            max_build_list.append([tech,            # TECHNOLOGY_TYPE
                                   tech_params[0],  # METHOD: 'PCT'/'ABS'
                                   tech_params[1],  # VALUE
                                   tech_params[2],  # START_YEAR
                                   tech_params[3]]) # END_YEAR
            '''
            max_build_temp = pd.DataFrame(columns=['TYPE', 
                                                   'METHOD', 
                                                   'MAX_BUILD',
                                                   'YEAR'])
            max_build_temp['YEAR'] = list(range(tech_params[2],
                                                tech_params[3]))
            max_build_temp['TYPE'] = tech
            max_build_temp['METHOD'] = tech_params[0]
            max_build_temp['MAX_BUILD'] = tech_params[1]
            max_build_df = pd.concat([max_build_df, max_build_temp],
                                     ignore_index=True)

        # Create a list of powerplant technologies
        tech_set = pd.read_csv(os.path.join(output_data_dir, 'TECHNOLOGY.csv'))
        pwr_tech_list = [x for x in list(tech_set['VALUE'])
                         if x.startswith('PWR')
                         ]
        
        # Create scaffold dataframe with all powerplant technologies for all years
        df_techs = pd.DataFrame(list(itertools.product(pwr_tech_list,
                                                       years)
                                     ),
                                columns = ['TECHNOLOGY', 
                                           'YEAR']
                                )
        
        # Filter out technologies for which a max. capacity investment has already
        # been set
        df_max_cap_inv = pd.read_csv(os.path.join(output_data_dir,
                                                  'TotalAnnualMaxCapacityInvestment.csv'))
        max_cap_inv_techs = list(df_max_cap_inv['TECHNOLOGY'].unique())
        df_techs = df_techs[~(df_techs['TECHNOLOGY'].isin(max_cap_inv_techs))]
        
        # Create dataframe of max capacity by technology
        if os.path.isfile(os.path.join(output_data_dir,
                                       'TotalAnnualMaxCapacity.csv')):
            df_max_cap = pd.read_csv(os.path.join(output_data_dir,
                                                  'TotalAnnualMaxCapacity.csv'))
            df_max_cap = df_max_cap.loc[df_max_cap['TECHNOLOGY'].str.startswith('PWR')]
        else:
            df_max_cap = pd.DataFrame(columns=['REGION',
                                               'TECHNOLOGY',
                                               'YEAR',
                                               'VALUE'])
        df_max_cap = pd.merge(left=df_techs, 
                              right=df_max_cap,
                              on=['TECHNOLOGY', 'YEAR'],
                              how='left')
        df_max_cap['REGION'] = region
        df_max_cap['TYPE'] = df_max_cap['TECHNOLOGY'].str[3:6]
        df_max_cap = pd.merge(left=df_max_cap, 
                              right=max_build_df,
                              on=['TYPE', 'YEAR'],
                              how='left')
        df_max_cap.loc[df_max_cap['METHOD'].isin(['ABS']),
                                  'VALUE'] = df_max_cap['MAX_BUILD']
        df_max_cap.dropna(inplace=True)
        df_max_cap.loc[df_max_cap['METHOD'].isin(['PCT']),
                                  'VALUE'] = df_max_cap['VALUE']*df_max_cap['MAX_BUILD']/100

        df_max_cap = df_max_cap[['REGION',
                                 'TECHNOLOGY',
                                 'YEAR',
                                 'VALUE']]
        df_max_cap = pd.concat([df_max_cap_inv,
                                df_max_cap],
                                ignore_index=True)
        df_max_cap.to_csv(os.path.join(output_data_dir, 
                                       "TotalAnnualMaxCapacityInvestment.csv"),
                                       index = None)


if __name__ == '__main__':
    main()
    logging.info(f'Max capacity limits sucessfully set')