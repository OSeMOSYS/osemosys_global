#!/usr/bin/env python
# coding: utf-8

# OSeMOSYS-PLEXOS global model: Powerplant data

# Import modules
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
from datetime import datetime
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import itertools
import urllib
import os
import yaml
from OPG_configuration import ConfigFile, ConfigPaths
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
import os
from scipy import spatial

def main():
    
    # CONFIGURATION PARAMETERS

    config_paths = ConfigPaths()
    config = ConfigFile('config')
    region_name = config.region_name

    input_data_dir = config_paths.input_data_dir
    output_data_dir = config_paths.output_data_dir

    custom_nodes_dir = config_paths.custom_nodes_dir

    # Check for custom nodes directory
    try:
        os.makedirs(custom_nodes_dir)
    except FileExistsError:
        pass

    cross_border_trade = config.get('crossborderTrade')
    model_start_year = config.get('startYear')
    model_end_year = config.get('endYear')
    years = list(range(model_start_year, model_end_year + 1))

    region_name = config.region_name
    tech_capacity = config.get('user_defined_capacity')
    custom_nodes = config.get('nodes_to_add')

    # Create output directory 
    if not os.path.exists(output_data_dir):
        os.makedirs(output_data_dir)

    #Checks whether PLEXOS-World 2015 data needs to be retrieved from the PLEXOS-World Harvard Dataverse.
    try:
        Open = open(os.path.join(input_data_dir,
                                 "PLEXOS_World_2015_Gold_V1.1.xlsx"))

    except IOError:
        urllib.request.urlretrieve("https://dataverse.harvard.edu/api/access/datafile/4008393?format=original&gbrecs=true" , 
                                   os.path.join(input_data_dir,
                                                "PLEXOS_World_2015_Gold_V1.1.xlsx")
                                   )

        Open = open(os.path.join(input_data_dir,
                                 "PLEXOS_World_2015_Gold_V1.1.xlsx")
                    )

    finally:
        Open.close()

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
                                           "weo_2020_powerplant_costs.csv")
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
    df_af = pd.read_csv(os.path.join(input_data_dir,
                                     "availability_factors.csv")
                        )
    
    if custom_nodes:
        df_custom_res_cap = pd.read_csv(os.path.join(input_data_dir,
                                                     "custom_nodes",
                                                     "residual_capacity.csv")
                                 )
    

    emissions = []

    # Technologies that will have 00 and 01 suffixes to represent PLEXOS 
    # historical values and future values 
    duplicate_techs = ['CCG', 'OCG', 'COA']

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
    if custom_nodes:
        node_list = list(df_gen_agg_node['node'].unique()) + custom_nodes
    else:
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
    tech_list = list(df_tech_code['code'])
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
    
    # Change IAR for URN value taken from PLEXOS to 2.2 (45%)
    df_eff.loc[df_eff['tech_code']=='URN', 'efficiency'] = 0.45

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

    if custom_nodes:
        df_res_cap_custom, custom_techs = custom_nodes_csv(custom_nodes, 
                                                           df_custom_res_cap, 
                                                           region_name, 
                                                           years, 
                                                           tech_list)
        df_res_cap = df_res_cap.append(df_res_cap_custom)

    # df_res_cap.to_csv(r"osemosys_global_model/data/ResidualCapacity.csv", index=None)
    df_res_cap.drop_duplicates(subset=['REGION','TECHNOLOGY','YEAR'],
                               keep='last',
                               inplace=True)
    df_res_cap.to_csv(os.path.join(output_data_dir, 
                                   "ResidualCapacity.csv"),
                      index=None)

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
            
    gen_locations = pd.merge(gen_locations, df_gen_2[['node', 'country_code', 'node_code', 'powerplant']], 
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
    if custom_nodes:
        node_list = list(df_gen_2['node_code'].unique()) + custom_nodes
    else:
        node_list = list(df_gen_2['node_code'].unique())
    # Add extra nodes which are not present in 2015 but will be by 2050
    for each_node in nodes_extra_list:
        if len(each_node) <= 6:
            node_list.append("".join(each_node.split('-')[1:]) + 'XX')
        else:
            node_list.append("".join(each_node.split('-')[1:]))

    master_fuel_list = list(df_gen_2['tech_code'].unique())
    master_fuel_list.append('CCS')

    mode_list = [1,2]

    df_ratios = pd.DataFrame(list(itertools.product(node_list,
                                                    master_fuel_list,
                                                    mode_list,
                                                    years)
                                  ),
                             columns = ['node_code', 'tech_code', 'MODE_OF_OPERATION', 'YEAR']
                             )

    df_ratios = createPwrTechs(df_ratios, duplicate_techs)

    thermal_fuel_list = ['COA',
                         'COG',
                         'OCG',
                         'CCG',
                         'PET',
                         'URN',
                         'OIL',
                         'OTH',
                         'CCS'
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
    
    # CCS Mode 1: Domestic COA
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['CCS'])),
               'FUEL'] = 'COA'+df_iar['TECHNOLOGY'].str[6:9]

    # CCS Mode 2: International COA
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['CCS'])),
               'FUEL'] = 'COAINT'

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
    
    df_iar.loc[df_iar['TECHNOLOGY'].str.startswith('PWRCCS'),
               'tech_average_iar'] = 3
    
    df_iar.loc[df_iar['VALUE'].isna(),
               'VALUE'] = df_iar['tech_average_iar']
    df_iar.drop_duplicates(inplace=True)

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

    # Create tech column from To and From Codes:
    df_trn_efficiencies = format_transmission_name(df_trn_efficiencies)

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
                                 'VALUE']]

    # Add iar for techs not using PLEXOS values 
    df_iar_newTechs = duplicatePlexosTechs(df_iar_final, duplicate_techs)
    for duplicate_tech in duplicate_techs:
        df_new_iar = newIar(df_iar_newTechs, duplicate_tech)
        df_iar_final = df_iar_final.append(df_new_iar)

    # Add oar for techs not using PLEXOS values 
    df_oar_newTechs = duplicatePlexosTechs(df_oar_final, duplicate_techs)
    df_oar_final = df_oar_final.append(df_oar_newTechs, ignore_index=True)

    #df_oar_final.to_csv(r"osemosys_global_model/data/OutputActivityRatio.csv", index = None)
    df_oar_final.drop_duplicates(inplace=True)
    df_oar_final.to_csv(os.path.join(output_data_dir,
                                     "OutputActivityRatio.csv"),
                        index=None)
    # df_iar_final.to_csv(r"osemosys_global_model/data/InputActivityRatio.csv", index = None)
    df_iar_final.to_csv(os.path.join(output_data_dir,
                                     "InputActivityRatio.csv"),
                        index=None)

    # ### Costs: Capital, fixed, and variable

    df_costs = pd.melt(df_weo_data, 
                       id_vars = ['technology', 'weo_region', 'parameter'], 
                       value_vars = ['2019', '2030', '2040'], 
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
                  'Other':'OTH',
                  'IGCC + CCS':'CCS',
                  'Coal* + CCS':'CCS',} # Added OIL, OTH, PET, WOF to WEO 2018

    df_costs = df_costs.loc[df_costs['technology'].isin(costs_dict.keys())]
    df_costs['technology_code'] = df_costs['technology'].replace(costs_dict)

    weo_regions_dict = dict([(k, v) 
                             for k, v 
                             in zip(df_weo_regions['technology_code'], 
                                    df_weo_regions['weo_region']
                                   )
                            ]
                           )
    
    # Get transmission costs
    df_trans_capex, df_trans_fix = get_transmission_costs()

    df_trans_capex['REGION'] = region_name
    df_trans_capex['YEAR'] = [years] * len(df_trans_capex)
    df_trans_capex = df_trans_capex.explode('YEAR')

    df_trans_fix['REGION'] = region_name
    df_trans_fix['YEAR'] = [years] * len(df_trans_fix)
    df_trans_fix = df_trans_fix.explode('YEAR')

    # Filter out techs that don't have activity ratios 
    df_trans_capex = df_trans_capex.loc[
        df_trans_capex['TECHNOLOGY'].isin(df_oar_final['TECHNOLOGY'])]
    df_trans_fix = df_trans_fix.loc[
        df_trans_fix['TECHNOLOGY'].isin(df_oar_final['TECHNOLOGY'])]

    # Create formatted CSVs
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
            df_costs_final = df_costs_final.merge(df_trans_capex, how='outer')
            df_costs_final.to_csv(os.path.join(output_data_dir, 
                                               "CapitalCost.csv"),
                                  index = None)
        if each_cost in ['O&M']:
            df_costs_final = df_costs_final.merge(df_trans_fix, how='outer')
            df_costs_final.to_csv(os.path.join(output_data_dir, 
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
    df_capact_final.to_csv(os.path.join(output_data_dir,
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
    else:
        df_crossborder_final = pd.DataFrame(columns=['REGION', 
                                                    'TECHNOLOGY',
                                                    'VALUE'])
    df_crossborder_final.to_csv(os.path.join(output_data_dir,
                                            "TotalTechnologyModelPeriodActivityUpperLimit.csv"),
                                index = None)


    # Create Operational Life data
    tech_code_dict_reverse = dict((v,k) for k,v in tech_code_dict.items())
    op_life_techs = list(set(list(df_iar_final['TECHNOLOGY'].unique()) + \
        list(df_oar_final['TECHNOLOGY'].unique())))
    op_life_techs = [tech for tech in op_life_techs if 
        (tech[0:3] == 'PWR') and (len(tech) == 13)]
    op_life_Out = []
    for op_life_tech in op_life_techs:
        op_life_tech_name = tech_code_dict_reverse[op_life_tech[3:6]]
        op_life_Out.append([
            region_name, 
            op_life_tech,
            op_life_dict[op_life_tech_name]])
    df_op_life_Out = pd.DataFrame(op_life_Out, columns = ['REGION', 'TECHNOLOGY', 'VALUE'])

    df_op_life_Out.to_csv(os.path.join(output_data_dir,
                                                "OperationalLife.csv"),
                                    index = None)

    # Create totalAnnualMaxCapacityInvestment data 

    # Do not allow capacity investment for all PWRxxxxxxxx00 technolgoies 
    max_cap_invest_techs = list(set(
        df_iar_final.loc[df_iar_final['TECHNOLOGY'].str.endswith('00')]['TECHNOLOGY'].tolist()))
    max_cap_invest_data = []
    for tech in max_cap_invest_techs:
        for year in years: 
            max_cap_invest_data.append([region_name, tech, year, 0])

    # Do not allow investment for all xxxABCxxxxxxx technologies
    no_investment_techs = config.get('no_invest_technologies')
    if not no_investment_techs:
        no_investment_techs = [] # Change from None type to empty list
    max_cap_invest_techs = list(set(df_iar_final.loc[
        df_iar_final['TECHNOLOGY'].str[3:6].isin(no_investment_techs)][
        'TECHNOLOGY'].tolist()))
    for tech in max_cap_invest_techs:
        for year in years: 
            max_cap_invest_data.append([region_name, tech, year, 0])
    
    # Save totalAnnualMaxCapacityInvestment
    df_max_cap_invest = pd.DataFrame(max_cap_invest_data,
                                    columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )       
    df_max_cap_invest.to_csv(os.path.join(output_data_dir, 
                                            'TotalAnnualMaxCapacityInvestment.csv'),
                                        index = None)
    
    df_min_cap_invest = pd.DataFrame(columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )
    df_min_cap_invest.to_csv(os.path.join(output_data_dir, 
                                            'TotalAnnualMinCapacityInvestment.csv'),
                                        index = None)

    # ## Create sets for TECHNOLOGIES, FUELS
    if custom_nodes:
        create_sets('TECHNOLOGY', df_oar_final, output_data_dir, custom_techs)
    else:
        create_sets('TECHNOLOGY', df_oar_final, output_data_dir, [])
    create_sets('FUEL', df_oar_final, output_data_dir, [])                             

    # ## Create set for YEAR, REGION, MODE_OF_OPERATION

    years_df = pd.DataFrame(years, columns = ['VALUE'])
    years_df.to_csv(os.path.join(output_data_dir, 
                                 "YEAR.csv"),
                    index = None)

    mode_list_df = pd.DataFrame(mode_list, columns = ['VALUE'])
    mode_list_df.to_csv(os.path.join(output_data_dir, 
                                     "MODE_OF_OPERATION.csv"),
                        index = None)

    regions_df = pd.DataFrame(columns = ['VALUE'])
    regions_df.loc[0] = region_name
    regions_df.to_csv(os.path.join(output_data_dir, 
                                   "REGION.csv"),
                      index = None)

    user_defined_capacity(region_name, years, output_data_dir, tech_capacity, op_life_dict)
    availability_factor(region_name, years, output_data_dir, df_af)


def create_sets(x, df, output_dir, custom_node_elements):
    """Creates a formatted otoole set csv 
    
    Arguments: 
        x = Name of set given by otoole as a string
        df = dataframe in extract set information from 
        output_dir = directory to write file to
    
    Returns: 
        None. Writes out csv file  
    
    Example:
        create_sets('TECHNOLOGY', df_oar_final, output_dir)
    """
    set_elements = list(df[x].unique()) + list(df[x].unique()) + list(custom_node_elements)
    set_elements = list(set(set_elements))
    set_elements = [x for x in set_elements if x != 'nan']
    set_elements.sort()
    set_elements_df = pd.DataFrame(set_elements, columns = ['VALUE'])
    return set_elements_df.to_csv(os.path.join(output_dir,
                                               str(x) + '.csv'
                                              ),
                                  index = None
                                 )

def duplicatePlexosTechs(df_in, techs):
    """Creates new technologies to replace PLEXOS technolgoies.
    
    New technologies will end in '01', while historical ones end in '00'
    
    Arguments: 
        df_in = dataframe in otoole and og formatting with a TECHNOLOGY column
        techs = List of technology triads to duplicate [CCG, HYD, ...]
    
    Returns: 
        df_out = dataframe with same columns as df_in. All tech names that include
        techs values will be returned with updated naming. Remaining technologies 
        are deleted 
    
    Example:
        df_out = duplicatePlexosTechs(df_in, ['CCG', 'OCG'])
        df_in['TECHNOLOGY'] = [PWRCCGAFGXX01, PWROCGAFGXX01, PWRHYDAFGXX01]
        df_out['TECHNOLOGY'] = [PWRCCGAFGXX02, PWROCGAFGXX02]
    """
    df_out = df_in.copy()
    # df_out = df_out.loc[df_out['TECHNOLOGY'].str[3:6].isin(techs)]
    df_out = df_out.loc[(df_out['TECHNOLOGY'].str[3:6].isin(techs)) & 
                        ~(df_out['TECHNOLOGY'].str.startswith('MIN'))]
    df_out['TECHNOLOGY'] = df_out['TECHNOLOGY'].str.slice_replace(start=11,
                                                                  stop=13,
                                                                  repl='01')
    return df_out

def createPwrTechs(df_in, techs):
    """Formats power generation technology name
    
    Adds a 'TECHNOLOGY' column to a dataframe with formatted power 
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
        df_in['tech_code'] = ('CCG', SPV, 'OCG', 'HYD')
        df_in['node_code'] = ('AGOXX', AGOXX, 'INDNP', 'INDNP')
        df_out = createPwrTechs(df_in, ['CCG', 'OCG'])
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

def newIar(df_in, tech):
    """Replaces the input activity ratio value with a hardcoded value 

    Arguments: 
        df = dataframe with a 'TECHNOLOGY' and 'VALUE' column
        tech = technology to replace iar for (CCG, HYD, SPV...)
    
    Returns: 
        df_out = same dataframe as df_in with a new values in 'VALUE'
    
    Example:
        df_out = newIar(df_in, 'CCG')
        df_out['TECHNOLOGY'] = [PWRCCGINDNP01, PWRCCGINDNW01]
        df_out['VALUE'] = [2, 2]
    """

    df_out = df_in.loc[df_in['TECHNOLOGY'].str[3:6] == tech]
    if tech == 'CCG':
        iar = 0.5
    elif tech == 'OCG':
        iar = 0.35
    elif tech == 'COA':
        iar = 0.33
    else: 
        logging.warning(f'Default IAR used for new {tech} power plants')
        iar = 1
    df_out['VALUE'] = round(1/iar, 3)
    return df_out

def get_transmission_costs():
    '''Gets electrical transmission capital and fixed cost per technology. 

    Both the capital costs and fixed cost are written out to avoid having 
    to read in the excel file twice 
    
    Returns: 
        df_capex: DataFrame with the columns 'TECHNOLOGY' and 'VALUE' 
            representing CAPITAL cost in millions of dollars per year. 
        df_fix: DataFrame with the columns 'TECHNOLOGY' and 'VALUE' 
            representing FIXED cost in millions of dollars per year. 
    '''

    # CONFIGURATION PARAMETERS

    config_paths = ConfigPaths()
    input_data_dir = config_paths.input_data_dir

    # Read in raw data
    df = pd.read_excel(os.path.join(input_data_dir,
                                                 "Costs Line expansion.xlsx"),
                                    sheet_name = 'Lines')

    # Drop unneeded columns
    df = df.drop(
        [
            "Line",
            "KM distance",
            "HVAC/HVDC/Subsea",
            "Losses",
            "Unnamed: 8",
            "Line Max Size (MW)",
            "Unnamed: 10",
            "Unnamed: 11",
            "Unnamed: 12",
            "Subsea lines",
            "Unnamed: 14"
        ],
        axis=1,
    )

    # Use to/from codes to create a TECHNOLOGY columns
    df = format_transmission_name(df)

    # Changes units
    # Raw given in dollars -> model units in million dollars
    df = df.rename(columns={'Annual FO&M (3.5% of CAPEX) ($2010 in $000)':'O&M'})
    df = df.rename(columns={'Build Cost ($2010 in $000)':'CAPEX'})
    df['O&M'] = df['O&M'].astype(float) / 1000
    df['CAPEX'] = df['CAPEX'].astype(float) / 1000
    df = df.round({'O&M': 3, 'CAPEX': 3, })

    # Separate out fixed and capex and return values 
    df_capex = df.drop(['O&M'], axis=1)
    df_capex = df_capex.rename(columns={'CAPEX':'VALUE'})
    df_fix = df.drop(['CAPEX'], axis=1)
    df_fix = df_fix.rename(columns={'O&M':'VALUE'})

    return df_capex, df_fix

def format_transmission_name(df):
    '''Formats PLEXOS transmission names into OSeMOSYS Global names.

    Args:
        :param df: Pandas DataFrame with columns 'From' and 'To' describing the 
               transmission from and to contries. ie. 
    
    Returns: 
        :param df: Same as df_in, except the 'From' and 'To' columns are replaced 
            with a single 'TECHNOLOGY' column holding OSeMOSYS Global 
            naming conventions 

    Example:
        df = pd.DataFrame(
            [[AF-COD, AF-COG, 0.001],
            [EU-AUT, EU-SVK, 0.004],
            [AS-LBN, AS-SYR, 0.006]], 
            columns = ['From', 'To', 'Losses']
        )
        pd.DataFrame(
            [[0.001,TRNCODXXCOGXX],
            [0.004,TRNAUTXXSVKXX],
            [0.006,TRNLBNXXSYRXX]] 
            columns = ['Losses','TECHNOLOGY'])'''

    # If from column has length 6 then it's the last three chars plus XX
    df.loc[df["From"].str.len() == 6, "From"] = (df["From"].str[3:6] + "XX")

    # If from column has length 9 then it's the 3:6 and 7:9 three chars plus XX
    df.loc[df["From"].str.len() == 9, "From"] = (
        df["From"].str[3:6] + df["From"].str[7:9])

    # If to column has length 6 then it's the last three chars plus XX
    df.loc[df["To"].str.len() == 6, "To"] = (df["To"].str[3:6] + "XX")

    # If to column has length 9 then it's the 3:6 and 7:9 three chars plus XX
    df.loc[df["To"].str.len() == 9, "To"] = (
        df["To"].str[3:6] + df["To"].str[7:9])

    # Combine From and To columns.
    df["TECHNOLOGY"] = ("TRN" + df["From"] + df["To"])

    # Drop to and from columns
    df = df.drop(["From", "To"], axis=1)

    return df

def user_defined_capacity(region, years, output_data_dir, tech_capacity, op_life_dict):
    """User-defined capacities are used when a specific technology must be 
    invested, for a given year and capacity. This is applied hrough the 
    parameter 'TotalAnnualMinCapacityInvestment'. 

    Args:
        region: From config file (e.g. 'GLOBAL')
        years: From config file (e.g. [2020, 2021, ... 2050])
        output_data_dir: Output directory set in config file
        tech_capacity: User-defined capacity in config file 
                       (e.g. TRNAGOXXCODXX: [5, 2030])

    Returns
        None
    """
    techCapacity = []
    tech_capacity_dict = {}
    
    if not tech_capacity is None:

        for tech, tech_params in tech_capacity.items():
            techCapacity.append([tech, tech_params[0], tech_params[1]])
            tech_capacity_dict[tech] = tech_params[2]
        tech_capacity_df = pd.DataFrame(techCapacity,
                                        columns=['TECHNOLOGY', 'VALUE', 'YEAR'])
        tech_capacity_df['REGION'] = region
        tech_capacity_df = tech_capacity_df[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]

        tech_set = pd.read_csv(os.path.join(output_data_dir, 'TECHNOLOGY.csv'))

        for each_tech in list(tech_capacity_df['TECHNOLOGY'].unique()):
            if each_tech not in list(tech_set['VALUE']):
        #        tech_capacity_df = tech_capacity_df.loc[~(tech_capacity_df['TECHNOLOGY'].isin([each_tech]))]
                tech_set = tech_set.append(pd.DataFrame({'VALUE':[each_tech]}))

        df_min_cap_inv = pd.read_csv(os.path.join(output_data_dir,
                                                  'TotalAnnualMinCapacityInvestment.csv'))
        df_min_cap_inv = df_min_cap_inv.append(tech_capacity_df)
        df_min_cap_inv.drop_duplicates(inplace=True)

        df_max_cap_inv = pd.read_csv(os.path.join(output_data_dir,
                                                  'TotalAnnualMaxCapacityInvestment.csv'))

        max_cap_techs = []
        for index, row in tech_capacity_df.iterrows():
            for each_year in years:
                if row['YEAR'] == each_year:
                    value = row['VALUE']
                if row['YEAR'] > each_year:
                    value = 0
                else:
                    if tech_capacity_dict[row['TECHNOLOGY']] in ['open']:
                        value = np.nan
                    elif tech_capacity_dict[row['TECHNOLOGY']] in ['fixed']:
                        value = 0
                max_cap_techs.append([row['REGION'],
                                      row['TECHNOLOGY'],
                                      each_year,
                                      value])
        max_cap_techs_df = pd.DataFrame(max_cap_techs,
                                        columns=['REGION',
                                                'TECHNOLOGY',
                                                'YEAR',
                                                'VALUE'])
        max_cap_techs_df.dropna(inplace=True)
        df_max_cap_inv = df_max_cap_inv.append(max_cap_techs_df)
        df_max_cap_inv.drop_duplicates(inplace=True)

        df_max_cap_inv.to_csv(os.path.join(output_data_dir,
                                        "TotalAnnualMaxCapacityInvestment.csv"),
                            index=None)
        # For technologies with start year before model start year, add to 
        # ResidualCapacity
        df_res_cap_ud = df_min_cap_inv.loc[df_min_cap_inv['YEAR'] < min(years)]
        df_res_cap_ud.rename(columns={'YEAR':'START_YEAR'},
                             inplace=True)
        df_res_cap_ud_final = pd.DataFrame(list(itertools.product(df_res_cap_ud['TECHNOLOGY'].unique(),
                                                                  years)
                                                ),
                                           columns = ['TECHNOLOGY',
                                                      'YEAR']
                                           )
        df_res_cap_ud_final = pd.merge(df_res_cap_ud_final,
                                       df_res_cap_ud,
                                       how='left',
                                       on=['TECHNOLOGY'])
        df_res_cap_ud_final['TECH'] = df_res_cap_ud_final['TECHNOLOGY'].str[3:6]
        df_res_cap_ud_final.loc[df_res_cap_ud_final['TECHNOLOGY'].str.contains('TRN'),
                                'TECH'] = 'TRN'
        df_res_cap_ud_final['OP_LIFE'] = df_res_cap_ud_final['TECH'].map(op_life_dict)
        df_res_cap_ud_final['END_YEAR'] = (df_res_cap_ud_final['OP_LIFE'] 
                                           + df_res_cap_ud_final['START_YEAR'])
        df_res_cap_ud_final = df_res_cap_ud_final.loc[df_res_cap_ud_final['YEAR'] 
                                                      >= df_res_cap_ud_final['START_YEAR']]
        df_res_cap_ud_final = df_res_cap_ud_final.loc[df_res_cap_ud_final['YEAR'] 
                                                      <= df_res_cap_ud_final['END_YEAR']]
        df_res_cap_ud_final['REGION'] = region
        df_res_cap_ud_final = df_res_cap_ud_final[['REGION',
                                                   'TECHNOLOGY',
                                                   'YEAR',
                                                   'VALUE']]
        df_res_cap = pd.read_csv(os.path.join(output_data_dir,
                                              'ResidualCapacity.csv'))
        df_res_cap = pd.concat([df_res_cap, df_res_cap_ud_final])
        df_res_cap.to_csv(os.path.join(output_data_dir,
                                       'ResidualCapacity.csv'),
                          index=None)
                
        # For technologies with start year at or after model start year, add to 
        # TotalAnnualMinCapacityInvestment      
        df_min_cap_inv = df_min_cap_inv.loc[df_min_cap_inv['YEAR'] >= min(years)]
        df_min_cap_inv.to_csv(os.path.join(output_data_dir,
                                        "TotalAnnualMinCapacityInvestment.csv"),
                            index=None)
        tech_set.to_csv(os.path.join(output_data_dir,
                                     "TECHNOLOGY.csv"),
                        index=None)
        
        # Add IAR and OAR for custom technologies
        df_iar = pd.read_csv(os.path.join(output_data_dir,
                                                  'InputActivityRatio.csv'))
        df_oar = pd.read_csv(os.path.join(output_data_dir,
                                                  'OutputActivityRatio.csv'))
        tech_list = list(tech_capacity_df['TECHNOLOGY'].unique())
        df_iar_custom = pd.DataFrame(list(itertools.product(tech_list,
                                                        [1, 2],
                                                        years)
                                      ),
                                 columns = ['TECHNOLOGY',
                                            'MODE_OF_OPERATION',
                                            'YEAR']
                                 )
        df_oar_custom = pd.DataFrame(list(itertools.product(tech_list,
                                                        [1, 2],
                                                        years)
                                      ),
                                 columns = ['TECHNOLOGY',
                                            'MODE_OF_OPERATION',
                                            'YEAR']
                                 )  
        # IAR in modes 1 and 2 are primary electricity commodity ('ELC*01') in 
        # node_from and node_to, respectively. 
        # OAR is the inverse of the above
        df_iar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==1,
                                        'FUEL'] = ('ELC' + 
                                                   df_iar_custom['TECHNOLOGY'].str[3:8] + 
                                                   '01')
        df_iar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==2,
                                        'FUEL'] = ('ELC' + 
                                                   df_iar_custom['TECHNOLOGY'].str[8:13] + 
                                                   '01')
        df_oar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==1,
                                        'FUEL'] = ('ELC' + 
                                                   df_oar_custom['TECHNOLOGY'].str[8:13] + 
                                                   '01')
        df_oar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==2,
                                        'FUEL'] = ('ELC' + 
                                                   df_oar_custom['TECHNOLOGY'].str[3:8] + 
                                                   '01')
        df_iar_custom['VALUE'] = 1
        df_oar_custom['VALUE'] = 0.9
        df_iar_custom['REGION'] = region
        df_oar_custom['REGION'] = region

        df_iar_custom = df_iar_custom[['REGION', 
                                       'TECHNOLOGY',
                                       'FUEL', 
                                       'MODE_OF_OPERATION',
                                       'YEAR', 
                                       'VALUE',]]
        df_oar_custom = df_oar_custom[['REGION', 
                                       'TECHNOLOGY',
                                       'FUEL', 
                                       'MODE_OF_OPERATION',
                                       'YEAR', 
                                       'VALUE',]]
        
        df_iar = pd.concat([df_iar, df_iar_custom])
        df_oar = pd.concat([df_oar, df_oar_custom])
        
        df_iar.drop_duplicates(inplace=True)
        df_iar.to_csv(os.path.join(output_data_dir,
                                   'InputActivityRatio.csv'),
                      index=None)
        df_oar.drop_duplicates(inplace=True)
        df_oar.to_csv(os.path.join(output_data_dir,
                                   'OutputActivityRatio.csv'),
                      index=None)
        # Add new fuels to FUEL set, if not already present
        fuel_set = pd.read_csv(os.path.join(output_data_dir, 'FUEL.csv'))
        fuel_list = []
        fuel_list = list(df_iar_custom['FUEL'].unique()) + list(df_oar_custom['FUEL'].unique())
        fuel_list = list(set(fuel_list))
        for each_fuel in fuel_list:
            if each_fuel not in list(fuel_set['VALUE']):
                fuel_set = fuel_set.append(pd.DataFrame({'VALUE':[each_fuel]}))

        fuel_set.to_csv(os.path.join(output_data_dir,
                                     "FUEL.csv"),
                        index=None)

        op_life = pd.read_csv(os.path.join(output_data_dir,
                                           'OperationalLife.csv'))
        op_life_custom = pd.DataFrame({'TECHNOLOGY': tech_list})

        op_life_custom.loc[op_life_custom['TECHNOLOGY'].str.contains('TRN'),
                           'VALUE'] = 60
        op_life_custom['REGION'] = region
        op_life_custom = op_life_custom[['REGION',
                                         'TECHNOLOGY',
                                         'VALUE']]

        op_life = pd.concat([op_life, op_life_custom])
        op_life.to_csv(os.path.join(output_data_dir,
                                    'OperationalLife.csv'),
                       index=None)
        # Add CapacityToActivityUnit for custom technologies
        cap_act = pd.read_csv(os.path.join(output_data_dir,
                                           'CapacityToActivityUnit.csv'))
        cap_act_custom = pd.DataFrame({'TECHNOLOGY': tech_list})
        cap_act_custom.loc[cap_act_custom['TECHNOLOGY'].str.contains('TRN'),
                           'VALUE'] = 31.536
        cap_act_custom.loc[cap_act_custom['TECHNOLOGY'].str.contains('PWR'),
                           'VALUE'] = 31.536
        cap_act_custom['REGION'] = region
        cap_act_custom = cap_act_custom[['REGION',
                                         'TECHNOLOGY',
                                         'VALUE']]
        cap_act = pd.concat([cap_act, cap_act_custom])
        cap_act.to_csv(os.path.join(output_data_dir,
                                    'CapacityToActivityUnit.csv'),
                       index=None)
        
        # Add CapitalCost for custom technologies
        cap_cost = pd.read_csv(os.path.join(output_data_dir,
                                            'CapitalCost.csv'))
        tech_list = list(tech_capacity_df['TECHNOLOGY'].unique())
        cap_cost_trn = pd.DataFrame(list(itertools.product(tech_list,
                                                           years)),
                                    columns = ['TECHNOLOGY',
                                               'YEAR'])
        cap_cost_trn.loc[cap_cost_trn['TECHNOLOGY'].str.contains('TRN'),
                         'VALUE'] = 1100
        cap_cost_trn.loc[cap_cost_trn['TECHNOLOGY'].str.contains('PWRTRN'),
                         'VALUE'] = 800
        cap_cost_trn['REGION'] = region
        cap_cost_trn = cap_cost_trn[['REGION',
                                     'TECHNOLOGY',
                                     'YEAR',
                                     'VALUE']]
        cap_cost = pd.concat([cap_cost, cap_cost_trn])
        cap_cost.to_csv(os.path.join(output_data_dir,
                                     'CapitalCost.csv'),
                        index=None)

def custom_nodes_csv(custom_nodes, df_custom, region, years, tech_list):
    '''Add custom nodes to the model for each relevant input parameter data csv.

    Args:
        df : Pandas DataFrame with columns 'From' and 'To' describing the 
                transmission from and to contries. ie. 
    
    Returns: 
        df_out : 

    '''
    df_param = pd.DataFrame(list(itertools.product(custom_nodes,
                                                   tech_list,
                                                   years)
                                  ),
                             columns = ['CUSTOM_NODE',
                                        'FUEL_TYPE',
                                        'YEAR']
                             )
    df_param['REGION'] = region
    df_custom = df_custom.groupby(['CUSTOM_NODE',
                                   'FUEL_TYPE',
                                   'START_YEAR',
                                   'END_YEAR'],
                                  as_index=False)['CAPACITY'].sum()
    df_param = pd.merge(df_param,
                        df_custom,
                        how='left',
                        on=['CUSTOM_NODE',
                            'FUEL_TYPE'])
    df_param['TECHNOLOGY'] = ('PWR' +
                              df_param['FUEL_TYPE'] + 
                              df_param['CUSTOM_NODE'] +
                              '01')
    technologies = df_param['TECHNOLOGY'].unique()
    df_param.dropna(inplace=True)
    df_param.drop_duplicates(inplace=True)
    df_param = df_param.loc[df_param['YEAR'] >= df_param['START_YEAR']]
    df_param = df_param.loc[df_param['YEAR'] <= df_param['END_YEAR']]
    df_param['VALUE'] = df_param['CAPACITY'].div(1000)
    df_param['REGION'] = region
    df_param = df_param[['REGION','TECHNOLOGY','YEAR','VALUE']]
    df_param = df_param.groupby(['REGION',
                                 'TECHNOLOGY',
                                 'YEAR'],
                                 as_index=False)['VALUE'].sum()

    return df_param, technologies


def availability_factor(region, 
                        years,
                        output_data_dir,
                        availability):
    
    af_dict = dict(zip(list(availability['technology']),
                       list(availability['value'])))
    
    df_tech = pd.read_csv(os.path.join(output_data_dir,
                                       'TECHNOLOGY.csv'))
    tech_list = [x for x in df_tech['VALUE']
                 if x.startswith('PWR')]
    df_af_final = pd.DataFrame(list(itertools.product(tech_list,
                                                      years)
                                    ),
                               columns = ['TECHNOLOGY', 'YEAR']
                               )
    df_af_final['TECH'] = df_af_final['TECHNOLOGY'].str[3:6]
    df_af_final['VALUE'] = df_af_final['TECH'].map(af_dict)
    df_af_final.dropna(inplace=True)
    df_af_final['REGION'] = region
    
    df_af_final = df_af_final[['REGION',
                               'TECHNOLOGY',
                               'YEAR',
                               'VALUE']]
    df_af_final.to_csv(os.path.join(output_data_dir,
                                    'AvailabilityFactor.csv'),
                       index=None)

if __name__ == "__main__":
    main()
    logging.info('Powerplant Data Created')
