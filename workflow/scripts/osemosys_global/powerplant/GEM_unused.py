"""Function that calculates existing, planned and to be retired capacities based on GEM for coal and gas. CURRENTLY UNUSED"""

import pandas as pd
import os
from scipy import spatial
import numpy as np

from constants import (
    start_year,
    input_data_dir
)

"""The output of the gem_cap function is currently unused. Implementation is being discussed so kept as standalone
function for the time being."""
def gem_cap(df_gen_base, tech_code_dict, op_life_dict):

    # ### Calculate planned capacities based on the Global Energy Observatory datasets
        
    # Add spatial mapping to link Global Energy Observatory region naming conventions with OSeMOSYS Global
    df_gem_regions = pd.read_csv(os.path.join(input_data_dir,
                                                  "gem_region_mapping.csv"), encoding = "ISO-8859-1")
    
    # pull locations from existing powerplants in PLEXOS-World dataset
    gen_locationsinput = pd.read_excel(os.path.join(input_data_dir, 
                                         "PLEXOS_World_2015_Gold_V1.1.xlsx"), 
                            sheet_name = "Attributes")
    
    gen_lat = gen_locationsinput.loc[(gen_locationsinput['class'].isin(['Battery', 'Generator'])) & 
                                                (gen_locationsinput['attribute'] == 'Latitude')].set_index('name')
                                                
    gen_long = gen_locationsinput.loc[(gen_locationsinput['class'].isin(['Battery', 'Generator'])) & 
                                                (gen_locationsinput['attribute'] == 'Longitude')].set_index('name')                                  
    
    gen_locations = pd.DataFrame(gen_long.index).merge(
        gen_lat['value'], left_on = 'name', right_index = True).merge(
            gen_long['value'], left_on = 'name', right_index = True)
            
    gen_locations = pd.merge(gen_locations, df_gen_base[['node', 'country_code', 'node_code', 'powerplant']], 
                             left_on = 'name', right_on = 'powerplant', how = 'inner')
    
    subcountries = gen_locations.loc[~gen_locations['node_code'].str.contains('XX')]['country_code'].str[:3].unique()
    
    # Set column name dictionaries for different Global Energy Monitor (gem) input datasets
    gem_coal_col = {'Country' : 'Country', 'Capacity (MW)' : 'VALUE', 
                    'Status' : 'Status', 'Year' : 'Year_built',
                    'RETIRED' : 'Year_retired', 'Planned Retire' : 'Year_retired_planned', 
                    'Latitude' : 'Latitude', 'Longitude' : 'Longitude'}
    
    gem_gas_col = {'Country' : 'Country', 'Capacity elec. (MW)' : 'VALUE', 
                    'Status' : 'Status', 'Start year' : 'Year_built',
                    'Retired year' : 'Year_retired', 'Planned retire' : 'Year_retired_planned', 
                    'Latitude' : 'Latitude', 'Longitude' : 'Longitude', 
                    'Technology' : 'Technology'}
    
    # Set technology dictionary to match with OSeMOSYS global technologies
    gem_gas_dict = {'CC' : 'CCG', 
                    'GT' : 'OCG',
                    'ICCC' : 'CCG',
                    'ISCC' : 'CCG',
                    'ST' : 'OCG',
                    'AFC' : 'CCG'}
    
    # Set which status criteria are used to filter datasets
    new_criteria = ['operating', 'proposed', 'announced', 'pre-permit', 'permitted', 'construction']# 'shelved' & cancelled' not included
    old_criteria = ['mothballed', 'retired', 'operating']# operating added because currently operating plants can already have an intended retirement year added
    
    # Import gem Datasets
    gem_coal = pd.read_excel(os.path.join(input_data_dir, 
                                                     'Global-Coal-Plant-Tracker-Jan-2022.xlsx'),
                                        sheet_name = 'Units', usecols = gem_coal_col.keys())
    
    
    gem_gas = pd.read_excel(os.path.join(input_data_dir, 
                                                     'Global-Gas-Plant-Tracker-Feb-2022.xlsx'),
                                        sheet_name = 'Gas Units', usecols = gem_gas_col.keys())
    
    # Add Technology columns and drop entries with no capacity values
    gem_coal.rename(columns = gem_coal_col, inplace = True)
    gem_coal['Technology'] = 'COA'
    gem_coal = gem_coal[pd.to_numeric(gem_coal['VALUE'], errors='coerce').notnull()]
    
    gem_gas.rename(columns = gem_gas_col, inplace = True)
    gem_gas['Technology'] = gem_gas['Technology'].map(gem_gas_dict)
    gem_gas = gem_gas[pd.to_numeric(gem_gas['VALUE'], errors='coerce').notnull()]
    
    # For entries in the gem dataset with no specified gas technology same assumptions are applied as with OSeMOSYS Global
    gem_gas.loc[
       (gem_gas["Technology"].isna()) & (gem_gas["VALUE"].astype(float) > 130),
       "Technology",
    ] = "CCG"
    gem_gas.loc[
        (gem_gas["Technology"].isna()) & (gem_gas["VALUE"].astype(float) <= 130),
        "Technology",
    ] = "OCG"
    
    # Combine different datasets
    gem_all = pd.concat([gem_coal, gem_gas])
    
    # Add spatial mapping
    gem_all = gem_all.merge(df_gem_regions, left_on = 'Country', right_on = 'gem_region')

    gem_concat = pd.DataFrame(columns = gem_all.columns)
    
    # Matches the lat/longs of plants in the gem datasets with lat/longs of the nearest plants in PLEXOS-World. 
    # The associated sub-country node of the nearest plant is assumed to be the node for the gem dataset entry. 
    for a in subcountries:
        
        gen_locations_sc = gen_locations.loc[gen_locations['country_code'].str.contains(a)].reset_index(drop = True)
        gem_all_sc = gem_all.loc[gem_all['country_code'] == a].reset_index(drop = True)
    
        source = spatial.KDTree(gen_locations_sc[['value_x', 'value_y']].to_numpy())
        
        output = pd.DataFrame(source.query([gem_all_sc[['Latitude', 'Longitude']
                                                                       ].to_numpy()])[1].transpose())
        
        gem_all_sc = gem_all_sc.merge(output, left_index = True, right_index = True).set_index(0)
        gem_all_sc = gem_all_sc.merge(gen_locations_sc[['node_code']], left_index = True, 
                                      right_index = True)
        
        gem_concat = pd.concat([gem_concat, gem_all_sc])
    
    # Adds matched sub-country entries to original df and sets node codes.    
    gem_all = gem_all.loc[~gem_all['country_code'].isin(subcountries)]
    gem_all = pd.concat([gem_all, gem_concat], axis = 0).reset_index(drop = True)
    gem_all['node_code'].fillna(gem_all['country_code'] + 'XX', inplace = True)
    
    # Filters datframe for new plants by year entry
    gem_all_new = gem_all.loc[(gem_all['Status'].isin(new_criteria)) & (gem_all['Year_built'].notna())
                              ].reset_index(drop = True)
    
    # Strips year entry to single value (last year taken if range is given e.g. 2025-)
    gem_all_new['YEAR'] = np.where(gem_all_new['Year_built'].astype(str).str.len() == 4, gem_all_new['Year_built'], 
                                   gem_all_new['Year_built'].str[-4:])
    
    # Drops non-existing or non-numerical entries (e.g. 'no year found') and entries from < model_start_year
    gem_all_new['YEAR'] = gem_all_new['YEAR'].apply(pd.to_numeric, errors = 'coerce')
    gem_all_new = gem_all_new[(gem_all_new['YEAR'].notna()) & 
                              (gem_all_new['YEAR'] > start_year)].reset_index()
    
    gem_all_retired = gem_all.loc[(gem_all['Status'].isin(old_criteria)) & (gem_all['Year_retired'].notna()) | 
                                  (gem_all['Year_retired_planned'].notna())]
    
    # Pulls retirement OR planned retirement year
    gem_all_retired['YEAR'] = np.where(gem_all_retired['Year_retired'].notna(), gem_all_retired['Year_retired'], 
                                       gem_all_retired['Year_retired_planned'])
    
    # Strips year entry to single value (last year taken if range is given e.g. 2025-2030)
    gem_all_retired['YEAR'] = np.where(gem_all_retired['YEAR'].astype(str).str.len().isin({4,6}), gem_all_retired['YEAR'], 
                                       gem_all_retired['YEAR'].str[-4:])
    
    # Drops non-existing or non-numerical entries (e.g. 'no year found') and entries from < model_start_year
    gem_all_retired['YEAR'] = gem_all_retired['YEAR'].apply(pd.to_numeric, errors = 'coerce')
    gem_all_retired = gem_all_retired[(gem_all_retired ['YEAR'].notna()) & 
                                      (gem_all_retired['YEAR'] > start_year)].reset_index()
    
    # Group values by technology, node & year
    gem_all_new_agg = gem_all_new.groupby(['node_code', 'Technology', 'YEAR'], 
                                              as_index=False)['VALUE'].sum()
    
    # Adds lifetime to planned capacities, calculates future retirement year and adds to retirement dataframe. 
    tech_code_dict_inv = {v: k for k, v in tech_code_dict.items()}
    gem_all_new_agg_oplife = gem_all_new_agg.copy()
    gem_all_new_agg_oplife['operational_life'] =  gem_all_new_agg_oplife['Technology'
                                                                         ].map(tech_code_dict_inv).map(op_life_dict)
    
    gem_all_new_agg_oplife['YEAR'] = gem_all_new_agg_oplife['YEAR'] + gem_all_new_agg_oplife['operational_life']
    
    gem_all_retired = pd.concat([gem_all_retired, gem_all_new_agg_oplife])
    
    gem_all_retired_agg = gem_all_retired.groupby(['node_code', 'Technology', 'YEAR'], 
                                              as_index=False)['VALUE'].sum()
