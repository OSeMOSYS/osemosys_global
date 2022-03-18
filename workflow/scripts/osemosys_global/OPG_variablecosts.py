#!/usr/bin/env python
# coding: utf-8

# # OSeMOSYS-PLEXOS global model: Powerplant data

# ### Import modules

# In[260]:

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import itertools
import os
import yaml
from OPG_configuration import ConfigFile, ConfigPaths
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# CONFIGURATION PARAMETERS

config_paths = ConfigPaths()
config = ConfigFile('config')

input_dir = config_paths.input_dir
input_data_dir = config_paths.input_data_dir
output_dir = config_paths.output_dir
output_data_dir = config_paths.output_data_dir


# ### Import data files and user input              

## Data for variable costs of fuels taken from World Bank Commodity Market Outlooks:
##    https://www.worldbank.org/en/research/commodity-markets
## Download the 'Charts and Data' zip file and extract the forecasts file (CMO-October-2021-forecasts.xlsx).
## Adjust filename below for forecasts done at a different date.

# Read in World Bank Commodity Price Outlook - we only want rows 87 - 91
# using 85 as the headers (years) and skipping the energy header...

#### REPLACE FROM BELOW ONCE FIGURED OUT...df_prices = pd.read_excel(
df_prices = pd.read_excel(
    os.path.join(input_data_dir,
    "CMO-October-2021-forecasts.xlsx"), header=1, skiprows=83, nrows=6
)

# Read in Technologies
df_techs = pd.read_csv(os.path.join(output_data_dir,'TECHNOLOGY.csv'))
df_trn_techs = df_techs.copy()

model_start_year = config.get('startYear')
model_end_year = config.get('endYear')
years = [range(model_start_year, model_end_year + 1)]

region_name = config.get('region')
emissions = []


# ### Filter technologies to keep only fuel production technologies (MIN)

df_techs = df_techs[df_techs.VALUE.str.contains('MIN')]

df_techs['REGION'] = region_name
df_techs['MODE_OF_OPERATION'] = 1
df_techs_temp = df_techs.copy()
df_techs_temp['MODE_OF_OPERATION'] = 2
df_techs = pd.concat([df_techs,df_techs_temp])
df_techs.rename(columns={'VALUE':'TECHNOLOGY'}, inplace=True)


df_techs['YEAR'] = [range(model_start_year, model_end_year+1)] * len(df_techs)

df_techs = df_techs.explode('YEAR')

df_techs.reset_index(drop=True, inplace=True)

# print(df_techs)


# ### Cleanup prices from CMO data

df_prices = df_prices.drop([0], axis=0)
df_prices = df_prices.drop(["Commodity"], axis=1)
columns = df_prices.columns.values
columns[0] = "YEAR"
df_prices = df_prices.loc[:, ~df_prices.columns.str.contains('^Unnamed', na=False)]

df_prices = df_prices.transpose()

# use commodity key as header
new_header = df_prices.iloc[0] #grab the first row for the header
df_prices = df_prices[1:] #take the data less the header row
new_header[1] = 'MINCOA'  # This was KFCOAL_AUS  KFCRUDE_PETRO  KFNGAS_EUR
new_header[2] = 'MINOIL'  # This was KFCRUDE_PETRO
new_header[3] = 'MINGAS'  # This was KFNGAS_EUR
df_prices.columns = new_header #set the header row as the df header

# drop units as we're doing $/PJ
df_prices = df_prices.drop("Unit", axis = 0)

# And convert to $/PJ
# Values taken from kylesconverter.com:
#  1 MT coal contains 29.31 PJ (1,000,000 tonnes coal) Original was $/mt ($/tonne I think).  We want $mill/PJ so divide by 29.31
#  1 MMBtu = 0.000001055056 PJ Original was MMBtu and we want $mill/pJ, so divide by 1.055056.
#  1 bbl = 0.00000612 PJ (Barrels of Oil) Original was $/bbl so divide by 0.00000612 and multiply by 1000000
# NOT ENTIRELY SURE I'M CALCULATING THE UNITS CORRECTLY (NOT SURE WHAT ORIGINAL UNITS WERE AS IT'S NOT STATED ANYWHERE)

# ORIGINAL UNITS:
#       KFCOAL_AUS KFCRUDE_PETRO KFNGAS_EUR KFNGAS_US KFNGAS_JP
#Unit       $/mt         $/bbl    $/mmbtu   $/mmbtu   $/mmbtu

df_prices['MINCOA'] = df_prices['MINCOA'] / 0.00002931
df_prices['MINOIL'] = df_prices['MINOIL'] / 6.12
df_prices['MINGAS'] = df_prices['MINGAS'] / 1.055056
df_prices['KFNGAS_US'] = df_prices['KFNGAS_US'] / 1.055056
df_prices['KFNGAS_JP'] = df_prices['KFNGAS_JP'] / 1.055056

# print(df_prices)

# Add in other fuels that are the same as those above:
df_prices['MINCOG'] = df_prices['MINGAS']  # Cogen is powered by gas
df_prices['MINOTH'] = df_prices['MINOIL']  # Other petroleum products are similar to oil
df_prices['MINPET'] = df_prices['MINOIL']  # Petroleum products are similar to oil

# And add in international prices that are 15% higher than the regular ones:
df_prices['INTCOA'] = df_prices['MINCOA'] * 0.15
df_prices['INTOIL'] = df_prices['MINOIL'] * 0.15
df_prices['INTGAS'] = df_prices['MINGAS'] * 0.15
df_prices['INTCOG'] = df_prices['MINCOG'] * 0.15
df_prices['INTOTH'] = df_prices['MINOTH'] * 0.15
df_prices['INTPET'] = df_prices['MINPET'] * 0.15

df_prices = df_prices.reindex(range(model_start_year, model_end_year+1))

df_prices = df_prices.apply(pd.to_numeric)

df_prices = df_prices.interpolate()

col = df_prices.columns
df_prices = pd.melt(df_prices.reset_index(), id_vars='index', value_vars=col)

df_prices = df_prices.rename(columns={'index': 'YEAR', 'YEAR': 'TEMPTECH'})

# print(df_prices)


# ### Map costs to technologies by region/country

## Need to create: REGION,TECHNOLOGY,MODE_OF_OPERATION,YEAR,VALUE

# Setup TEMPTECH column for merge
df_techs['TEMPTECH'] = df_techs['TECHNOLOGY'].str[0:6]

# International costs are identified by the INT in MINCOAINT, so make these INTCOA temporarily
df_techs.loc[df_techs.TECHNOLOGY.str[6:9]=='INT', 'TEMPTECH'] = 'INT'+df_techs['TECHNOLOGY'].str[3:6]

df_varcost = pd.merge(df_techs, df_prices, on=['YEAR', 'TEMPTECH'])

df_varcost = df_varcost.rename(columns={'value': 'VALUE'})

df_varcost = df_varcost.drop(['TEMPTECH'], axis=1)

# print(df_varcost)

## AND THEN DROP ANY ITEMS THAT HAVE A NAN AS THESE ARE DEFAULT (0) VALUES
## BUT THIS IS ALREADY DONE FOR US BY THE MERGE
#df_techs.dropna(subset = ["VALUE"], inplace=True)
#print(df_techs)



# Get transmission variable costs 

df_trn_techs = df_trn_techs.loc[
    (df_trn_techs['VALUE'].str.startswith('TRN')) &
    (df_trn_techs['VALUE'].str.len() == 13)] # keeps TRNxxxxxxxxxx techs

df_trn_varcosts = df_trn_techs.copy()
df_trn_varcosts = df_trn_varcosts.rename(columns={'VALUE':'TECHNOLOGY'})
df_trn_varcosts['REGION'] = region_name
df_trn_varcosts['YEAR'] = years * len(df_trn_varcosts)
df_trn_varcosts['MODE_OF_OPERATION'] = [[1,2] for x in range(len(df_trn_varcosts))]
df_trn_varcosts = df_trn_varcosts.explode('MODE_OF_OPERATION')
df_trn_varcosts = df_trn_varcosts.explode('YEAR')

# Hardcode in transmission variable cost parameter following PLEXOS of $4/MWh
# $4/MWh * 1Wh/3600J * 1000000000 MJ/1PJ * 1M$/$1000000 
trn_varcost = 4 / 3.6

df_trn_varcosts['VALUE'] = round(trn_varcost, 4)

# Merge with mining variable costs 
df_varcost = pd.concat([df_varcost, df_trn_varcosts])

# ### Write out variablecost.csv

df_varcosts_final = df_varcost[['REGION', 
                       'TECHNOLOGY', 
                       'MODE_OF_OPERATION',
                       'YEAR', 
                       'VALUE']]

df_varcosts_final.to_csv(os.path.join(output_data_dir,'VariableCost.csv'), mode='w', header=True, index = None)

logging.info('Variable Costs Completed')
