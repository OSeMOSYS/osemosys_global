# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 20:23:21 2022

@author: mbrinkerink
"""
 
# %%

def createPwrTechs(df_in, techs):
    """Adds a 'TECHNOLOGY' column to a dataframe with formatted power 
    generation technology (PWR) names. Names are formatted so the suffix 
    for plants in the 'techs' argument list will have 00, while everything 
    else will have an 01 suffix

    Arguments: 
    df_in = dataframe with a 'tech_codes' and 'node_codes' column
    tList = List of technology triads to have 00 suffix [CCG, HYD, ...]
    
    Returns: 
    df_out = same dataframe as df_in, except with a 'TECHNOLOGY' column added 
    to the end 
    
    Example:
    df_out = createPwrTechs(df_in, ['CCG', 'OCG'])
    df_in['tech_code'] = ('CCG', SPV, 'OCG', 'HYD')
    df_in['node_code'] = ('AGOXX', AGOXX, 'INDNP', 'INDNP')
    df_out['TECHNOLOGY'] = [PWRCCGAFGXX00, PWRSPVAFGXX01, PWROCGINDNP00, PWRHYDINDNP01]
    """
    df_out = df_in.copy()
    for t in techs:
        df_out.loc[df_out['tech_code'] == t, 'tech_suffix'] = '00'
    df_out['tech_suffix'] = df_out['tech_suffix'].fillna('01')
    df_out['TECHNOLOGY'] = ('PWR' + 
                            df_out['tech_code'] + 
                            df_out['node_code'] + 
                            df_out['tech_suffix']
                            )
    df_out = df_out.drop('tech_suffix', axis = 1)
    return df_out

# Import modules
import pandas as pd
from datetime import datetime
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import itertools
import urllib
import os
import yaml
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
import os

yaml_file = open(os.path.join(os.path.dirname(__file__), '../../..',
                              'config/config.yaml'))
parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

input_dir = os.path.join(os.path.dirname(__file__), '../../..',
    parsed_yaml_file.get('inputDir'))
input_data_dir = os.path.join(input_dir, 'data')

output_dir = os.path.join(os.path.dirname(__file__), '../../..', 
    parsed_yaml_file.get('outputDir'))
output_data_dir =  os.path.join(output_dir, 'data')

df = pd.read_excel(os.path.join(input_data_dir,
                                "PLEXOS_World_2015_Gold_V1.1.xlsx"), 
                   sheet_name = "Properties")

df_dict = pd.read_excel(os.path.join(input_data_dir, 
                                     "PLEXOS_World_2015_Gold_V1.1.xlsx"), 
                        sheet_name = "Memberships")

df_dict = df_dict[df_dict["parent_class"] == "Generator"].rename(
    {"parent_object": "powerplant"}, axis=1
    )
df_weo_data = pd.read_csv(os.path.join(input_data_dir,
                                       "weo_2018_powerplant_costs.csv")
                          )
df_op_life = pd.read_csv(os.path.join(input_data_dir,
                                      "operational_life.csv")
                         )
df_tech_code = pd.read_csv(os.path.join(input_data_dir,
                                        "naming_convention_tech.csv")
                           )
df_trn_efficiencies = pd.read_excel(os.path.join(input_data_dir,
                                                 "Costs Line expansion.xlsx"),
                                    sheet_name = 'Interface'
                                    )
df_weo_regions = pd.read_csv(os.path.join(input_data_dir,
                                          "weo_region_mapping.csv")
                             )

model_start_year = parsed_yaml_file.get('startYear')
model_end_year = parsed_yaml_file.get('endYear')
years = list(range(model_start_year, 
                   model_end_year + 1))

region_name = parsed_yaml_file.get('region')
emissions = []

# Technologies that will have 00 and 01 suffixes to represent PLEXOS 
# historical values and future values 
duplicate_techs = ['CCG', 'OCG']

# Create main generator table
gen_cols_1 = ["child_class", "child_object", "property", "value"]
df_gen = df[gen_cols_1]
df_gen = df_gen[df_gen["child_class"] == "Generator"]
df_gen.rename(columns={"child_object": "powerplant"}, inplace=True)
df_gen.drop("child_class", axis=1, inplace=True)
df_gen = pd.pivot_table(df_gen,
                        index="powerplant",
                        columns="property",
                        values="value",
                        aggfunc=np.sum,
                        fill_value=0,
                       )
df_gen["total_capacity"] = (df_gen["Max Capacity"].astype(float)) * (
    df_gen["Units"].astype(int)
)

gen_cols_2 = ["Commission Date", "Heat Rate", "Max Capacity", "total_capacity"]
df_gen_2 = df_gen[gen_cols_2]

## Compile dataframe with powerplants, nodes, and fuels
df_dict_fuel = df_dict[df_dict["collection"] == "Fuels"]
df_dict_fuel = df_dict_fuel[["powerplant", "child_object"]]
df_dict_nodes = df_dict[df_dict["collection"] == "Nodes"]
df_dict_nodes = df_dict_nodes[["powerplant", "child_object"]]
df_dict_2 = pd.merge(df_dict_fuel, df_dict_nodes, how="outer", on="powerplant")


## Merge original generator dataframe with nodes and fuels
df_gen_2 = pd.merge(df_gen_2, df_dict_2, how="outer", on="powerplant")
df_gen_2.rename(
    {"child_object_x": "fuel", "child_object_y": "node"}, axis=1, inplace=True
)

## Extract start year from Commission Date
df_gen_2["Commission Date"] = (pd.TimedeltaIndex(df_gen_2["Commission Date"].astype(int),
                                            unit='d') + 
                           datetime(1900, 1, 1))
df_gen_2["start_year"] = df_gen_2["Commission Date"].dt.year
df_gen_2.drop("Commission Date", axis=1, inplace=True)

## Calculate efficiency from heat rate. Units of heat rate in MJ/kWh
df_gen_2["efficiency"] = 3.6 / df_gen_2["Heat Rate"].astype(float)
df_gen_2.drop("Heat Rate", axis=1, inplace=True)

## Calcluate years of operation from start year until 2015
df_gen_2["years_of_operation"] = model_start_year - df_gen_2["start_year"]

## Fix blank spaces in 'fuels' columns. Appearing for 'Oil' powerplants in certain countries
df_gen_2.loc[df_gen_2["fuel"].isna(), "fuel"] = (
    df_gen_2["node"].str.split("-").str[:2].str.join("-")
    + " "
    + df_gen_2["powerplant"].str.split("_", expand=True)[1]
)


## Create column for technology
df_gen_2["technology"] = df_gen_2["powerplant"].str.split("_").str[1]
df_gen_2["technology"] = df_gen_2["technology"].str.title()


## Divide Gas into CCGT and OCGT based on max capacity
df_gen_2.loc[
    (df_gen_2["technology"] == "Gas") & (df_gen_2["Max Capacity"].astype(float) > 130),
    "technology",
] = "Gas-CCGT"
df_gen_2.loc[
    (df_gen_2["technology"] == "Gas") & (df_gen_2["Max Capacity"].astype(float) <= 130),
    "technology",
] = "Gas-OCGT"


# Create table with aggregated capacity  
df_gen_agg_node = df_gen_2[df_gen_2['start_year']<=model_start_year]
df_gen_agg_node = df_gen_agg_node.groupby(['node', 'technology'], 
                                          as_index=False)['total_capacity'].sum()
df_gen_agg_node = df_gen_agg_node.pivot(index='node', 
                                        columns='technology', 
                                        values='total_capacity').fillna(0).reset_index()

df_gen_agg_node.drop('Sto', axis=1, inplace=True) # Drop 'Sto' technology. Only for USA.

# Add extra nodes which exist in 2050 but are not in the 2015 data
node_list = list(df_gen_agg_node['node'].unique())
nodes_extra_df = pd.DataFrame(columns=['node'])
nodes_extra_list = ['AF-SOM',
                    'AF-TCD',
                    'AS-TLS',
                    'EU-MLT',
                    'NA-BLZ',
                    'NA-HTI',
                    'SA-BRA-J1',
                    'SA-BRA-J2',
                    'SA-BRA-J3',
                    'SA-SUR',]
nodes_extra_df['node'] = nodes_extra_list

df_gen_agg_node = df_gen_agg_node.append(nodes_extra_df,
                                         ignore_index=True,
                                         sort='False').fillna(0).sort_values(by='node').set_index('node').round(2)
#df_gen_agg_node.to_csv(r'output/test_output_2.csv')


# Add region and country code columns
df_gen_2['region_code'] = df_gen_2['node'].str[:2]
df_gen_2['country_code'] = df_gen_2['node'].str[3:]


# ### Add operational life column
op_life_dict = dict(zip(list(df_op_life['tech']),
                        list(df_op_life['years'])))

df_gen_2['operational_life'] = df_gen_2['technology'].map(op_life_dict)
df_gen_2['retirement_year_data'] = (df_gen_2['operational_life'] 
                                    + df_gen_2['start_year'])
df_gen_2['retirement_diff'] = ((df_gen_2['years_of_operation'] 
                               - df_gen_2['operational_life'])/
                               df_gen_2['operational_life'])

''' Set retirement year based on years of operation. 
If (years of operation - operational life) is more than 50% of 
operational life, set retirement year
'''
df_gen_2.loc[df_gen_2['retirement_diff'] >= 0.5, 
             'retirement_year_model'] = 2025
df_gen_2.loc[(df_gen_2['retirement_diff'] < 0.5) &
             (df_gen_2['retirement_diff'] > 0), 
             'retirement_year_model'] = 2030
df_gen_2.loc[df_gen_2['retirement_diff'] <= 0, 
             'retirement_year_model'] = df_gen_2['retirement_year_data']

#df_gen_2.to_csv(r'output/test_output_3.csv')


# ### Add naming convention
tech_code_dict = dict(zip(list(df_tech_code['tech']),
                          list(df_tech_code['code'])))
df_gen_2['tech_code'] = df_gen_2['technology'].map(tech_code_dict)

df_gen_2.loc[df_gen_2['node'].str.len() <= 6, 
             'node_code'] = (df_gen_2['node'].
                             str.split('-').
                             str[1:].
                             str.join("") +
                             'XX')
df_gen_2.loc[df_gen_2['node'].str.len() > 6, 
             'node_code'] = (df_gen_2['node'].
                             str.split('-').
                             str[1:].
                             str.join("")
                            )

df_gen_2 = df_gen_2.loc[~df_gen_2['tech_code'].isna()]


# ### Calculate average InputActivityRatio by node+technology and only by technology
df_eff = df_gen_2[['node_code',
                   'efficiency',
                   'tech_code']]

# Change IAR for CSP value taken from PLEXOS to 1.0
df_eff.loc[df_eff['tech_code']=='CSP', 'efficiency'] = 1

# Average efficiency by node and technology
df_eff_node = df_eff.groupby(['tech_code',
                              'node_code'],
                             as_index = False).agg('mean')

df_eff_node['node_average_iar'] = ((1 / df_eff_node['efficiency']).
                                   round(2))

df_eff_node.drop('efficiency', 
                 axis = 1, 
                 inplace = True)

# Average efficiency by technology
df_eff_tech = df_eff.groupby('tech_code',
                             as_index = False).agg('mean')

df_eff_tech['tech_average_iar'] = ((1 / df_eff_tech['efficiency']).
                                   round(2))

df_eff_tech.drop('efficiency', 
                 axis = 1, 
                 inplace = True)

# %%

# ### Calculate residual capacity
res_cap_cols = [
    "node_code",
    "tech_code",
    "total_capacity",
    "start_year",
    "retirement_year_model",
]

df_res_cap = df_gen_2[res_cap_cols]

for each_year in range(model_start_year, model_end_year+1):
    df_res_cap[str(each_year)] = 0

df_res_cap = pd.melt(
    df_res_cap,
    id_vars=res_cap_cols,
    value_vars=[x for x in df_res_cap.columns if x not in res_cap_cols],
    var_name="model_year",
    value_name="value",
)
df_res_cap["model_year"] = df_res_cap["model_year"].astype(int)
df_res_cap.loc[
    (df_res_cap["model_year"] >= df_res_cap["start_year"])
    & (df_res_cap["model_year"] <= df_res_cap["retirement_year_model"]),
    "value",
] = df_res_cap["total_capacity"]

df_res_cap = df_res_cap.groupby(
    ["node_code", "tech_code", "model_year"], as_index=False
)["value"].sum()

# Add column with naming convention
df_res_cap = createPwrTechs(df_res_cap, duplicate_techs)

# Convert total capacity from MW to GW
df_res_cap['value'] = df_res_cap['value'].div(1000)


df_res_cap_plot = df_res_cap[['node_code', 
                             'tech_code', 
                             'model_year', 
                             'value']]

# Rename 'model_year' to 'year' and 'total_capacity' to 'value'
df_res_cap.rename({'model_year': 'YEAR',
                   'value': 'VALUE'},
                  inplace=True,
                  axis=1)
# Drop 'tech_code' and 'node_code'
df_res_cap.drop(['tech_code', 'node_code'], 
                inplace=True, 
                axis=1)

# Add 'REGION' column and fill 'GLOBAL' throughout
df_res_cap['REGION'] = region_name

# Reorder columns
df_res_cap = df_res_cap[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]

# %% BELOW HERE IS NEW

from scipy import spatial

# Add spatial mapping
df_gem_regions = pd.read_csv(os.path.join(input_data_dir,
                                              "gem_region_mapping.csv"), encoding = "ISO-8859-1")

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
        
gen_locations = pd.merge(gen_locations, df_gen_2[['node', 'country_code', 'node_code', 'powerplant']], 
                         left_on = 'name', right_on = 'powerplant', how = 'inner')

subcountries = gen_locations.loc[~gen_locations['node_code'].str.contains('XX')]['country_code'].str[:3].unique()

# %% 

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

# %%
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
                          (gem_all_new['YEAR'] > model_start_year)].reset_index()

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
                                  (gem_all_retired['YEAR'] > model_start_year)].reset_index()

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