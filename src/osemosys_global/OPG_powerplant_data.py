#!/usr/bin/env python
# coding: utf-8

# OSeMOSYS-PLEXOS global model: Powerplant data

# Import modules
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import itertools
import urllib
import os
import yaml


#Read in information from YAML file
yaml_file = open("config.yaml")
parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

input_dir = parsed_yaml_file.get('inputDir')
output_dir = parsed_yaml_file.get('outputDir') + 'data/'

cross_border_trade = parsed_yaml_file.get('crossborderTrade')

#Checks whether PLEXOS-World 2015 data needs to be retrieved from the PLEXOS-World Harvard Dataverse.
try:
    #Open = open(r"../../data/PLEXOS_World_2015_Gold_V1.1.xlsx")
    Open = open(os.path.join(input_dir,
                             "PLEXOS_World_2015_Gold_V1.1.xlsx"))

except IOError:
    urllib.request.urlretrieve("https://dataverse.harvard.edu/api/access/datafile/4008393?format=original&gbrecs=true" , 
                               os.path.join(input_dir,
                                            "PLEXOS_World_2015_Gold_V1.1.xlsx")
                               )

    Open = open(os.path.join(input_dir,
                             "PLEXOS_World_2015_Gold_V1.1.xlsx")
                )

finally:
    Open.close()

df = pd.read_excel(os.path.join(input_dir,
                                "PLEXOS_World_2015_Gold_V1.1.xlsx"), 
                   sheet_name = "Properties")

df_dict = pd.read_excel(os.path.join(input_dir, 
                                     "PLEXOS_World_2015_Gold_V1.1.xlsx"), 
                        sheet_name = "Memberships")

df_dict = df_dict[df_dict["parent_class"] == "Generator"].rename(
    {"parent_object": "powerplant"}, axis=1
    )
df_weo_data = pd.read_csv(os.path.join(input_dir,
                                       "weo_2018_powerplant_costs.csv")
                          )
df_op_life = pd.read_csv(os.path.join(input_dir,
                                      "operational_life.csv")
                         )
df_tech_code = pd.read_csv(os.path.join(input_dir,
                                        "naming_convention_tech.csv")
                           )
df_trn_efficiencies = pd.read_excel(os.path.join(input_dir,
                                                 "Costs Line expansion.xlsx"),
                                    sheet_name = 'Interface'
                                    )
df_weo_regions = pd.read_csv(os.path.join(input_dir,
                                          "weo_region_mapping.csv")
                             )

model_start_year = parsed_yaml_file.get('startYear')
model_end_year = parsed_yaml_file.get('endYear')
years = list(range(model_start_year, 
                   model_end_year + 1))

region_name = parsed_yaml_file.get('region')
emissions = []

# Create 'output' directory if it doesn't exist 
import os
#if not os.path.exists('osemosys_global_model/data'):
#    os.makedirs('osemosys_global_model/data')

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
df_gen_2["Commission Date"] = pd.to_datetime(df_gen_2["Commission Date"])
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
df_res_cap['node_code'] = df_res_cap['node_code']
df_res_cap['tech'] = ('PWR' + 
                      df_res_cap['tech_code'] + 
                      df_res_cap['node_code'] + '01'
                     )
# Convert total capacity from MW to GW
df_res_cap['value'] = df_res_cap['value'].div(1000)


df_res_cap_plot = df_res_cap[['node_code', 
                             'tech_code', 
                             'model_year', 
                             'value']]

# Rename 'model_year' to 'year' and 'total_capacity' to 'value'
df_res_cap.rename({'tech': 'TECHNOLOGY',
                   'model_year': 'YEAR',
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

# df_res_cap.to_csv(r"osemosys_global_model/data/ResidualCapacity.csv", index=None)
df_res_cap.to_csv(os.path.join(output_dir, 
                               "ResidualCapacity.csv"),
                  index=None)

'''
# ### Interactive visualisation of residual capacity by node

import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes = True)
from ipywidgets import interact, interactive, fixed, interact_manual, Layout
import ipywidgets as widgets
#importing plotly and cufflinks in offline mode
import plotly as py
#import plotly.graph_objs as go
import cufflinks
import plotly.offline as pyo
from plotly.offline import plot, iplot, init_notebook_mode
pyo.init_notebook_mode()
cufflinks.go_offline()
cufflinks.set_config_file(world_readable=True, theme='white')

color_codes = pd.read_csv(r'data\color_codes.csv', encoding='latin-1')
color_dict = dict([(n,c) for n,c in zip(color_codes.tech, color_codes.colour)])

def f(node):
    df_plot = df_res_cap_plot.loc[df_res_cap_plot['node_code']==node]
    df_plot.drop('node_code', 
                     axis = 1, 
                     inplace = True)
    df_plot = df_plot.pivot_table(index='model_year',
                                  columns='tech_code',
                                  values='value',
                                  aggfunc='sum').reset_index()
    
    
    #plt.figure(figsize=(10, 10), dpi= 80, facecolor='w', edgecolor='k')
    #ax = sns.barplot(df_plot)
    return df_plot.iplot(x = 'model_year',
                         kind = 'bar', 
                         barmode = 'stack',
                         xTitle = 'Year',
                         yTitle = 'Gigawatts',
                         color=[color_dict[x] for x in df_plot.columns if x != 'model_year'],
                         title = 'Residual Capacity',
                         showlegend = True)

interact(f,
         node=widgets.Dropdown(options = (df_res_cap_plot['node_code']
                                          .unique()
                                         )
                              )
        )
'''

# ### Add input and output activity ratios

# Create master table for activity ratios 
node_list = list(df_gen_2['node_code'].unique())

# Add extra nodes which are not present in 2015 but will be by 2050
for each_node in nodes_extra_list:
    if len(each_node) <= 6:
        node_list.append("".join(each_node.split('-')[1:]) + 'XX')
    else:
        node_list.append("".join(each_node.split('-')[1:]))

master_fuel_list = list(df_gen_2['tech_code'].unique())

mode_list = [1,2]

df_ratios = pd.DataFrame(list(itertools.product(node_list,
                                                master_fuel_list,
                                                mode_list,
                                                years)
                              ),
                         columns = ['node_code', 'tech_code', 'MODE_OF_OPERATION', 'YEAR']
                         )

df_ratios['TECHNOLOGY'] = ('PWR' + 
                           df_ratios['tech_code'] + 
                           df_ratios['node_code'] + '01'
                          )

thermal_fuel_list = ['COA',
                     'COG',
                     'OCG',
                     'CCG',
                     'PET',
                     'URN',
                     'OIL',
                     'OTH'
                    ]

thermal_fuel_list_iar = ['COA',
                         'COG',
                         'PET',
                         'URN',
                         'OIL',
                         'OTH'
                        ]

renewables_list = ['BIO',
                   'GEO',
                   'HYD',
                   'SPV', 
                   'CSP',
                   'WAS',
                   'WAV',
                   'WON', 
                   'WOF']


# #### OutputActivityRatio - Power Generation Technologies
df_oar = df_ratios.copy()
mask = df_oar['TECHNOLOGY'].apply(lambda x: x[3:6] in thermal_fuel_list)
df_oar['FUEL'] = 0
df_oar['FUEL'][mask] = 1
df_oar = df_oar.loc[~((df_oar['MODE_OF_OPERATION'] > 1) &
                      (df_oar['FUEL'] == 0))]
df_oar['FUEL'] = ('ELC' + 
                  df_oar['TECHNOLOGY'].str[6:11] + 
                  '01'
                 )
df_oar['VALUE'] = 1

# Add 'REGION' column and fill 'GLOBAL' throughout
df_oar['REGION'] = region_name

# Select columns for final output table
df_oar_final = df_oar[['REGION', 
                 'TECHNOLOGY',
                 'FUEL',                  
                 'MODE_OF_OPERATION',
                 'YEAR', 
                 'VALUE',]]

# Don't write yet - we'll write the IAR and OAR at the end...
# df_oar_final.to_csv(r"output/OutputActivityRatio.csv", index = None)


# #### InputActivityRatio - Power Generation Technologies
# Copy OAR table with all columns to IAR
df_iar = df_oar.copy()

df_iar['FUEL'] = 0

# Deal with GAS techs first...  OCG and CCG
# OCG Mode 1: Domestic GAS
df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
           (df_iar['TECHNOLOGY'].str[3:6].isin(['OCG'])),
           'FUEL'] = 'GAS'+df_iar['TECHNOLOGY'].str[6:9]
# OCG Mode 2: International GAS
df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
           (df_iar['TECHNOLOGY'].str[3:6].isin(['OCG'])),
           'FUEL'] = 'GASINT'

# CCG Mode 1: Domestic GAS
df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
           (df_iar['TECHNOLOGY'].str[3:6].isin(['CCG'])),
           'FUEL'] = 'GAS'+df_iar['TECHNOLOGY'].str[6:9]

# CCG Mode 2: International GAS
df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
           (df_iar['TECHNOLOGY'].str[3:6].isin(['CCG'])),
           'FUEL'] = 'GASINT'

# For non-GAS thermal fuels, domestic fuel input by country in mode 1 and 
# 'international' fuel input in mode 2
df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
           (df_iar['TECHNOLOGY'].str[3:6].isin(thermal_fuel_list_iar)),
           'FUEL'] = df_iar['TECHNOLOGY'].str[3:9]

df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
           (df_iar['TECHNOLOGY'].str[3:6].isin(thermal_fuel_list_iar)),
           'FUEL'] = df_iar['TECHNOLOGY'].str[3:6] + 'INT'

# For renewable fuels, input by node in mode 1
df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
           (df_iar['TECHNOLOGY'].str[3:6].isin(renewables_list)),
           'FUEL'] = df_iar['TECHNOLOGY'].str[3:11]

# Remove mode 2 when not used
df_iar = df_iar.loc[df_iar['FUEL'] != 0]

# Join efficiency columns: one with node and technology average, and the
# other with technology average
df_iar = df_iar.join(df_eff_node.set_index(['tech_code', 'node_code']), 
                     on=['tech_code', 'node_code'])

df_iar = df_iar.join(df_eff_tech.set_index('tech_code'), 
                     on='tech_code')

# When available, choose node and technology average. Else, 
# choose technology average
df_iar['VALUE'] = df_iar['node_average_iar']
df_iar.loc[df_iar['VALUE'].isna(),
           'VALUE'] = df_iar['tech_average_iar']

# Add 'REGION' column and fill 'GLOBAL' throughout
df_iar['REGION'] = region_name

# Select columns for final output table
df_iar_final = df_iar[['REGION', 
                       'TECHNOLOGY',
                       'FUEL',  
                       'MODE_OF_OPERATION',
                       'YEAR', 
                       'VALUE',]]

# Don't write this yet - we'll write both IAR and OAR at the end...
# df_iar_final.to_csv(r"output/InputActivityRatio.csv", index = None)


# #### OutputActivityRatios - Upstream

thermal_fuels = ['COA',
                 'COG',
                 'GAS',
                 'PET',
                 'URN',
                 'OIL',
                 'OTH'
                    ]

# We have to create a technology to produce every fuel that is input into any of the power technologies:

df_oar_upstream = df_iar_final.copy()

# All mining and resource technologies have an OAR of 1...
df_oar_upstream['VALUE'] = 1

# Renewables - set the technology as RNW + FUEL
df_oar_upstream.loc[df_oar_upstream['FUEL'].str[0:3].isin(renewables_list),
           'TECHNOLOGY'] = 'RNW'+df_oar_upstream['FUEL']

# If the fuel is a thermal fuel, we need to create the OAR for the mining technology... BUT NOT FOR THE INT FUELS...
df_oar_upstream.loc[df_oar_upstream['FUEL'].str[0:3].isin(thermal_fuels) & ~(df_oar_upstream['FUEL'].str[3:6] == "INT"),
           'TECHNOLOGY'] = 'MIN'+df_oar_upstream['FUEL']

# Above should get all the outputs for the MIN technologies, but we need to adjust the mode 2 ones to just the fuel code (rather than MINCOAINT)
df_oar_upstream.loc[df_oar_upstream['MODE_OF_OPERATION']==2,
           'TECHNOLOGY'] = 'MIN'+df_oar_upstream['FUEL'].str[0:3]+df_oar_upstream['TECHNOLOGY'].str[6:9]
df_oar_upstream.loc[df_oar_upstream['MODE_OF_OPERATION']==2,
           'FUEL'] = df_oar_upstream['FUEL'].str[0:3]

# Now remove the duplicate fuels that the above created (because there's now a COA for each country, not each region, and GAS is repeated twice for each region as well):
df_oar_upstream.drop_duplicates(keep='first',inplace=True)

# Now we have to create the MINXXXINT technologies.  They are all based on the MODE_OF_OPERATION == 2:
df_oar_int = pd.DataFrame(df_oar_upstream.loc[df_oar_upstream['MODE_OF_OPERATION'] == 2, :])

# At this point we should have only the internationally traded fuels since they're all mode 2.  So we can make the tech MINXXXINT and that's that.
df_oar_int['TECHNOLOGY'] = 'MIN'+df_oar_int['FUEL']+'INT'
# And rename the fuel to XXXINT
df_oar_int['FUEL'] = df_oar_int['FUEL']+'INT'
df_oar_int['MODE_OF_OPERATION'] = 1  # This is probably not strictly necessary as long as they're always the same in and out...

# and de-duplicate this list:
df_oar_int.drop_duplicates(keep='first',inplace=True)


# #### Input Activity Ratios - Upstream

# All we need to do is take in the thermal fuels for the MINXXXINT technologies.  This already exists as df_oar_int with the XXINT fuel so we can simply copy that:
df_iar_int = df_oar_int.copy()
df_iar_int['FUEL'] = df_iar_int['FUEL'].str[0:3]


# #### Downstream Activity Ratios

# Build transmission system outputs

df_iar_trn = df_oar_final.copy()

# Change the technology name to PWRTRNXXXXX
df_iar_trn["TECHNOLOGY"] = "PWRTRN" + df_iar_trn["FUEL"].str[3:8]
# Make all modes of operation 1
df_iar_trn["MODE_OF_OPERATION"] = 1
# And remove all the duplicate entries
df_iar_trn.drop_duplicates(keep="first", inplace=True)

# OAR for transmission technologies is IAR, but the fuel is 02 instead of 01:
df_oar_trn = df_iar_trn.copy()
df_oar_trn["FUEL"] = df_oar_trn["FUEL"].str[0:8] + "02"

# Build international transmission system from original input data, but for Line rather than Generator:
int_trn_cols = ["child_class", "child_object", "property", "value"]
df_int_trn = df[int_trn_cols]
df_int_trn = df_int_trn[df_int_trn["child_class"] == "Line"]

# For IAR and OAR we can drop the value:
df_int_trn = df_int_trn.drop(["child_class", "value"], axis=1)

# Create MofO column based on property:
df_int_trn["MODE_OF_OPERATION"] = 1
df_int_trn.loc[df_int_trn["property"] == "Min Flow", "MODE_OF_OPERATION"] = 2

# Use the child_object column to build the technology names:
df_int_trn["codes"] = df_int_trn["child_object"].str.split(pat="-")

# If there are only two locations, then the node is XX
df_int_trn.loc[df_int_trn["codes"].str.len() == 2, "TECHNOLOGY"] = (
    "TRN" + df_int_trn["codes"].str[0] + "XX" + df_int_trn["codes"].str[1] + "XX"
)
# If there are four locations, the node is already included
df_int_trn.loc[df_int_trn["codes"].str.len() == 4, "TECHNOLOGY"] = (
    "TRN"
    + df_int_trn["codes"].str[0]
    + df_int_trn["codes"].str[1]
    + df_int_trn["codes"].str[2]
    + df_int_trn["codes"].str[3]
)
# If there are three items, and the last item is two characters, then the second item is an XX:
df_int_trn.loc[
    (df_int_trn["codes"].str.len() == 3) & (df_int_trn["codes"].str[2].str.len() == 2),
    "TECHNOLOGY",
] = (
    "TRN"
    + df_int_trn["codes"].str[0]
    + "XX"
    + df_int_trn["codes"].str[1]
    + df_int_trn["codes"].str[2]
)
# If there are three items, and the last item is three characters, then the last item is an XX:
df_int_trn.loc[
    (df_int_trn["codes"].str.len() == 3) & (df_int_trn["codes"].str[2].str.len() == 3),
    "TECHNOLOGY",
] = (
    "TRN"
    + df_int_trn["codes"].str[0]
    + df_int_trn["codes"].str[1]
    + df_int_trn["codes"].str[2]
    + "XX"
)

# Set the value (of either IAR or OAR) to 1
df_int_trn["VALUE"] = 1
df_int_trn["REGION"] = region_name

df_int_trn = df_int_trn.drop(["property", "child_object", "codes"], axis=1)
df_int_trn["YEAR"] = model_start_year
# Add in the years:
df_temp = df_int_trn.copy()
for year in range(model_start_year + 1, model_end_year + 1):
    df_temp["YEAR"] = year
    df_int_trn = df_int_trn.append(df_temp)

df_int_trn = df_int_trn.reset_index(drop=True)


# Now create the input and output activity ratios
df_int_trn_oar = df_int_trn.copy()
df_int_trn_iar = df_int_trn.copy()

# IAR Mode 1 is input from first country:
df_int_trn_iar.loc[df_int_trn_iar["MODE_OF_OPERATION"] == 1, "FUEL"] = (
    "ELC" + df_int_trn_iar["TECHNOLOGY"].str[3:8] + "02"
)
# IAR Mode 2 is input from second country:
df_int_trn_iar.loc[df_int_trn_iar["MODE_OF_OPERATION"] == 2, "FUEL"] = (
    "ELC" + df_int_trn_iar["TECHNOLOGY"].str[8:13] + "02"
)

# OAR Mode 2 is output to first country:
df_int_trn_oar.loc[df_int_trn_oar["MODE_OF_OPERATION"] == 2, "FUEL"] = (
    "ELC" + df_int_trn_oar["TECHNOLOGY"].str[3:8] + "01"
)
# OAR Mode 1 is out to the second country:
df_int_trn_oar.loc[df_int_trn_oar["MODE_OF_OPERATION"] == 1, "FUEL"] = (
    "ELC" + df_int_trn_oar["TECHNOLOGY"].str[8:13] + "01"
)

# Drop unneeded columns
df_trn_efficiencies = df_trn_efficiencies.drop(
    [
        "Interface",
        "KM distance",
        "HVAC/HVDC/Subsea",
        "Build Cost ($2010 in $000)",
        "Build cost + FOM (3.5% of Capex per year)",
        "Unnamed: 8",
        "Interface is per MW",
        "Unnamed: 10",
        "Unnamed: 11",
        "Unnamed: 12",
        "Subsea lines",
        "Unnamed: 14"
    ],
    axis=1,
)

# Drop NaN values
df_trn_efficiencies = df_trn_efficiencies.dropna(subset=["From"])

# Create To and From Codes:
# If from column has length 6 then it's the last three chars plus XX
df_trn_efficiencies.loc[df_trn_efficiencies["From"].str.len() == 6, "From"] = (
    df_trn_efficiencies["From"].str[3:6] + "XX"
)
# If from column has length 9 then it's the 3:6 and 7:9 three chars plus XX
df_trn_efficiencies.loc[df_trn_efficiencies["From"].str.len() == 9, "From"] = (
    df_trn_efficiencies["From"].str[3:6] + df_trn_efficiencies["From"].str[7:9]
)
# If from column has length 6 then it's the last three chars plus XX
df_trn_efficiencies.loc[df_trn_efficiencies["To"].str.len() == 6, "To"] = (
    df_trn_efficiencies["To"].str[3:6] + "XX"
)
# If from column has length 9 then it's the 3:6 and 7:9 three chars plus XX
df_trn_efficiencies.loc[df_trn_efficiencies["To"].str.len() == 9, "To"] = (
    df_trn_efficiencies["To"].str[3:6] + df_trn_efficiencies["To"].str[7:9]
)

# Combine From and To columns.
# If the From is earlier in the alphabet the technology is in order, add tech with mode 1.
df_trn_efficiencies["TECHNOLOGY"] = ("TRN" + df_trn_efficiencies["From"] + df_trn_efficiencies["To"])

# Drop to and from columns
df_trn_efficiencies = df_trn_efficiencies.drop(["From", "To"], axis=1)

# Rename column 'VALUES'
df_trn_efficiencies = df_trn_efficiencies.rename(columns={"Losses": "VALUE"})

# And adjust OAR values to be output amounts vs. losses:
df_trn_efficiencies['VALUE'] = 1.0 - df_trn_efficiencies['VALUE']

# and add values into OAR matrix
df_int_trn_oar = df_int_trn_oar.drop(["VALUE"], axis=1)
df_int_trn_oar = pd.merge(
    df_int_trn_oar, df_trn_efficiencies, how="outer", on="TECHNOLOGY"
)


# #### Output IAR and OAR

# In[ ]:


# Combine the pieces from above and output to csv:

df_oar_final = df_oar_final.append(df_oar_upstream) # add upstream production technologies
df_oar_final = df_oar_final.append(df_oar_int) # Add in path through international markets
df_oar_final = df_oar_final.append(df_oar_trn) # Add in domestic transmission
df_oar_final = df_oar_final.append(df_int_trn_oar) # Add in international transmission

# Select columns for final output table
df_oar_final = df_oar_final.dropna()
df_oar_final = df_oar_final[['REGION', 
                             'TECHNOLOGY',
                             'FUEL',  
                             'MODE_OF_OPERATION',
                             'YEAR', 
                             'VALUE',]]

df_iar_final = df_iar_final.append(df_iar_int) # Add in path through international markets
df_iar_final = df_iar_final.append(df_iar_trn) # Add in domestic transmission
df_iar_final = df_iar_final.append(df_int_trn_iar) # Add in international transmission

# Select columns for final output table
df_iar_final = df_iar_final.dropna()
df_iar_final = df_iar_final[['REGION', 
                             'TECHNOLOGY',
                             'FUEL',  
                             'MODE_OF_OPERATION',
                             'YEAR', 
                             'VALUE',]]

#df_oar_final.to_csv(r"osemosys_global_model/data/OutputActivityRatio.csv", index = None)
df_oar_final.to_csv(os.path.join(output_dir,
                                 "OutputActivityRatio.csv"),
                    index=None)
# df_iar_final.to_csv(r"osemosys_global_model/data/InputActivityRatio.csv", index = None)
df_iar_final.to_csv(os.path.join(output_dir,
                                 "InputActivityRatio.csv"),
                    index=None)


# ### Costs: Capital, fixed, and variable

df_costs = pd.melt(df_weo_data, 
                   id_vars = ['technology', 'weo_region', 'parameter'], 
                   value_vars = ['2017', '2030', '2040'], 
                   var_name = ['YEAR'])
df_costs['parameter'] = df_costs['parameter'].str.split('\r\n').str[0]
df_costs['value'] = df_costs['value'].replace({'n.a.':0})
df_costs['value'] = df_costs['value'].astype(float) 
df_costs = df_costs.pivot_table(index = ['technology', 'parameter', 'YEAR'], 
                                columns = 'weo_region', 
                                values = 'value').reset_index()
df_costs['AS_average'] = (df_costs['China'] + 
                            df_costs['India'] + 
                            df_costs['Japan'] + 
                            df_costs['Middle East']).div(4)
df_costs['NA_average'] = (df_costs['United States'])
df_costs['SA_average'] = (df_costs['Brazil'])
df_costs['Global_average'] = (df_costs['Africa'] +
                              df_costs['Brazil'] + 
                              df_costs['Europe'] +
                              df_costs['China'] + 
                              df_costs['India'] + 
                              df_costs['Japan'] + 
                              df_costs['Middle East'] +
                              df_costs['Russia'] +
                              df_costs['United States']).div(9)
df_costs = pd.melt(df_costs, 
                   id_vars = ['technology', 'parameter', 'YEAR'], 
                   value_vars = [x 
                                 for x 
                                 in df_costs.columns 
                                 if x not in ['technology', 'parameter', 'YEAR']
                                ]
                  )
df_costs['YEAR'] = df_costs['YEAR'].astype(int)
costs_dict = {'Biomass - waste incineration - CHP':'WAS',
              'Biomass Power plant':'BIO', 
              'CCGT':'CCG', 
              'CCGT - CHP':'COG', 
              'Concentrating solar power':'CSP',
              'Gas turbine':'OCG',
              'Geothermal':'GEO', 
              'Hydropower - large-scale':'HYD',
              'Marine':'WAV',
              'Nuclear':'URN', 
              'Solar photovoltaics - Large scale':'SPV', 
              'Steam Coal - SUBCRITICAL':'COA',
              'Steam Coal - SUPERCRITICAL':'COA', 
              'Steam Coal - ULTRASUPERCRITICAL':'COA',
              'Wind onshore':'WON',
              'Wind offshore':'WOF',
              'Petroleum':'PET',
              'Oil':'OIL',
              'Other':'OTH',} # Added OIL, OTH, PET, WOF to WEO 2018

df_costs = df_costs.loc[df_costs['technology'].isin(costs_dict.keys())]
df_costs['technology_code'] = df_costs['technology'].replace(costs_dict)

weo_regions_dict = dict([(k, v) 
                         for k, v 
                         in zip(df_weo_regions['technology_code'], 
                                df_weo_regions['weo_region']
                               )
                        ]
                       )

for each_cost in ['Capital', 'O&M']:
    df_costs_temp = df_costs.loc[df_costs['parameter'].str.contains(each_cost)]
    df_costs_temp.drop(['technology', 'parameter'], 
                       axis = 1, 
                       inplace = True)
    df_costs_final = df_oar_final[['REGION',
                                   'TECHNOLOGY',
                                   'YEAR'
                                  ]]
    df_costs_final['YEAR'] = df_costs_final['YEAR'].astype(int)
    df_costs_final = df_costs_final.drop_duplicates()
    df_costs_final = (df_costs_final
                      .loc[(df_costs_final['TECHNOLOGY']
                            .str.startswith('PWR')
                           ) & 
                           (~df_costs_final['TECHNOLOGY']
                            .str.contains('TRN')
                           )
                          ]
                     )
    df_costs_final['technology_code'] = df_costs_final['TECHNOLOGY'].str[3:6]
    df_costs_final['weo_region'] = df_costs_final['TECHNOLOGY'].str[6:9]
    df_costs_final['weo_region'] = (df_costs_final['weo_region']
                                         .replace(weo_regions_dict))
    
    df_costs_final = pd.merge(df_costs_final, 
                              df_costs_temp, 
                              on = ['technology_code', 'weo_region', 'YEAR'], 
                              how = 'left'
                             )
    df_costs_final.drop(['technology_code', 'weo_region'], 
                        axis = 1, 
                        inplace = True)
    df_costs_final = df_costs_final.fillna(-9)
    df_costs_final = pd.pivot_table(df_costs_final, 
                                    index = ['REGION', 'YEAR'], 
                                    columns = 'TECHNOLOGY', 
                                    values = 'value').reset_index()
    df_costs_final = df_costs_final.replace([-9],[np.nan])
    #df_costs_final.set_index(['REGION', 'YEAR'], 
    #                         inplace = True)
    
    
    df_costs_final = df_costs_final.interpolate(method = 'linear', 
                                                limit_direction='forward').round(2)
    df_costs_final = df_costs_final.interpolate(method = 'linear', 
                                                limit_direction='backward').round(2)
    df_costs_final = pd.melt(df_costs_final, 
                             id_vars = ['REGION', 'YEAR'], 
                             value_vars = [x for x in df_costs_final.columns
                                           if x not in ['REGION', 'YEAR']
                                          ],
                             var_name = 'TECHNOLOGY',
                             value_name = 'VALUE'
                            )
    df_costs_final = df_costs_final[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]
    df_costs_final = df_costs_final[~df_costs_final['VALUE'].isnull()]
    
    if each_cost in ['Capital']:
        df_costs_final.to_csv(os.path.join(output_dir, 
                                           "CapitalCost.csv"),
                              index = None)
    if each_cost in ['O&M']:
        df_costs_final.to_csv(os.path.join(output_dir, 
                                           "FixedCost.csv"),
                              index = None)


# Create CapacityToActivityUnit csv
df_capact_final = df_oar_final[['REGION',
                                'TECHNOLOGY'
                                ]]
df_capact_final = df_capact_final.drop_duplicates()
df_capact_final = (df_capact_final
                   .loc[(df_capact_final['TECHNOLOGY']
                         .str.startswith('PWR')
                         ) | 
                        (df_capact_final['TECHNOLOGY']
                         .str.contains('TRN')
                         )
                        ]
                   )
df_capact_final['VALUE'] = 31.536
df_capact_final.to_csv(os.path.join(output_dir,
                                    "CapacityToActivityUnit.csv"),
                       index = None)

# Set cross-border trade to 0 if False
if not cross_border_trade:
    df_crossborder_final = df_oar_final[['REGION',
                                        'TECHNOLOGY'
                                        ]]
    df_crossborder_final = df_crossborder_final.drop_duplicates()
    df_crossborder_final = (df_crossborder_final
                            .loc[df_crossborder_final['TECHNOLOGY']
                                .str.startswith('TRN')
                                ]
                            )
    df_crossborder_final['VALUE'] = 0
    df_crossborder_final.to_csv(os.path.join(output_dir,
                                            "TotalTechnologyModelPeriodActivityUpperLimit.csv"),
                                index = None)

# ## Create sets for TECHNOLOGIES, FUELS

def create_sets(x):
    set_elements = list(df_iar_final[x].unique()) + list(df_oar_final[x].unique())
    set_elements = list(set(set_elements))
    set_elements.sort()
    set_elements_df = pd.DataFrame(set_elements, columns = ['VALUE'])
    return set_elements_df.to_csv(os.path.join(output_dir,
                                               str(x) + '.csv'
                                              ),
                                  index = None
                                 )

create_sets('TECHNOLOGY')
create_sets('FUEL')                             


# ## Create set for YEAR, REGION, MODE_OF_OPERATION

years_df = pd.DataFrame(years, columns = ['VALUE'])
years_df.to_csv(os.path.join(output_dir, 
                             "YEAR.csv"),
                index = None)

mode_list_df = pd.DataFrame(mode_list, columns = ['VALUE'])
mode_list_df.to_csv(os.path.join(output_dir, 
                                 "MODE_OF_OPERATION.csv"),
                    index = None)

regions_df = pd.DataFrame(columns = ['VALUE'])
regions_df.loc[0] = region_name
regions_df.to_csv(os.path.join(output_dir, 
                               "REGION.csv"),
                  index = None)


# ## Create set for EMISSION 

emissions_df = pd.DataFrame(emissions, columns = ['VALUE'])
emissions_df.to_csv(os.path.join(output_dir, 
                                 "EMISSION.csv"),
                index = None)

