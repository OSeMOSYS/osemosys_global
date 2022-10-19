#!/usr/bin/env python
# coding: utf-8

# # OSeMOSYS-PLEXOS global model: TS-dependent parameters

# ### Import modules

# In[1]:


import pandas as pd
import datetime
import numpy as np
import itertools
import seaborn as sns; sns.set()
import matplotlib
import matplotlib.pyplot as plt
import urllib
import os
from OPG_configuration import ConfigFile, ConfigPaths
import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# ### Input data files and user input

# CONFIGURATION PARAMETERS

config_paths = ConfigPaths()
config = ConfigFile('config')

input_dir = config_paths.input_dir
input_data_dir = config_paths.input_data_dir
output_dir = config_paths.output_dir
output_data_dir = config_paths.output_data_dir

# helper functions

def apply_timeshift(x, timeshift):
    '''Applies timeshift to organize dayparts.
    
    Args:
        x = Value between 0-24
        timeshift = value offset from UTC (-11 -> +12)'''

    x += timeshift
    if x > 23:
        return x - 24
    elif x < 0:
        return x + 24
    else:
        return x

# Checks whether PLEXOS-World 2015 data needs to be retrieved from the PLEXOS-World Harvard Dataverse.
try:
    # Open = open(r'data/All_Demand_UTC_2015.csv')
    Open = open(os.path.join(input_data_dir,
                             'All_Demand_UTC_2015.csv')
                )    
    # demand_df = pd.read_csv(r'data/All_Demand_UTC_2015.csv' , encoding='latin-1')
    demand_df = pd.read_csv(os.path.join(input_data_dir,
                                         'All_Demand_UTC_2015.csv'),
                            encoding='latin-1')

except IOError:
    urllib.request.urlretrieve ('https://dataverse.harvard.edu/api/access/datafile/3985039?format=original&gbrecs=true', 
                                os.path.join(input_data_dir,
                                             'All_Demand_UTC_2015.csv')
                                )

    demand_df = pd.read_csv(os.path.join(input_data_dir,
                                         'All_Demand_UTC_2015.csv'),
                            encoding='latin-1')

seasons_raw = config.get('seasons')
seasonsData = []
for s, months in seasons_raw.items():
    for month in months:
        seasonsData.append([month, s]) 
seasons_df = pd.DataFrame(seasonsData, 
    columns = ['month', 'season'])
seasons_df = seasons_df.sort_values(by = ['month']).reset_index(drop = True)

dayparts_raw = config.get('dayparts')
daypartData = []
for dp, hr in dayparts_raw.items():
    daypartData.append([dp, hr[0], hr[1]])
dayparts_df = pd.DataFrame(daypartData, 
    columns = ['daypart', 'start_hour', 'end_hour'])
timeshift = config.get('timeshift')
dayparts_df['start_hour'] = dayparts_df['start_hour'].map(lambda x: apply_timeshift(x, timeshift))
dayparts_df['end_hour'] = dayparts_df['end_hour'].map(lambda x: apply_timeshift(x, timeshift))

daytype_included = config.get('daytype')
model_start_year = config.get('startYear')
model_end_year = config.get('endYear')
years = list(range(model_start_year, model_end_year+1))


# In[4]:


csp_df = pd.read_csv(os.path.join(input_data_dir,
                                  'CSP 2015.csv'),
                     encoding='latin-1')
csp_df.name = 'CSP'

spv_df = pd.read_csv(os.path.join(input_data_dir,
                                  'SolarPV 2015.csv'),
                     encoding='latin-1')
spv_df.name = 'SPV'

nodes = ['-'.join(x.split('-')[1:])
         for x in spv_df.columns
         if x
         not in ['Datetime']]
regions = [x
           for x in spv_df.columns
           if x
           not in ['Datetime']]

node_region_dict = dict(zip(nodes,
                            regions))

hyd_df = pd.read_csv(os.path.join(input_data_dir,
                                  'Hydro_Monthly_Profiles (15 year average).csv'),
                     encoding='latin-1')
hyd_df = hyd_df.loc[hyd_df['NAME'].str.endswith('Capacity Scaler')]
hyd_df['NAME'] = (hyd_df['NAME']
                  .str.split('_')
                  .str[0])
# Drop Brazil transmission nodes J1, J2, J3
brazil_j_nodes = ['BRA-J1', 'BRA-J2', 'BRA-J3']
hyd_df = hyd_df.loc[~hyd_df['NAME'].isin(brazil_j_nodes)]
hyd_df = hyd_df.set_index('NAME').T.reset_index()
hyd_df.rename(columns={'index': 'MONTH'},
              inplace=True)
hyd_df['MONTH'] = (hyd_df['MONTH']
                   .str.replace('M', '')
                   .astype(int))

hyd_df_processed = pd.DataFrame(columns=['Datetime'])
hyd_df_processed['Datetime'] = spv_df['Datetime']
hyd_df_processed['MONTH'] = (hyd_df_processed['Datetime']
                             .str.split('/')
                             .str[1]
                             .astype(int))
hyd_df_processed = pd.merge(hyd_df_processed,
                            hyd_df,
                            how='left',
                            on='MONTH')
hyd_df_processed.drop(columns='MONTH',
                      inplace=True)
hyd_df_processed.rename(columns=node_region_dict,
                        inplace=True)
hyd_df_processed.name = 'HYD'

won_df = pd.read_csv(os.path.join(input_data_dir,
                                  'Won 2015.csv'),
                     encoding='latin-1')
won_df.name = 'WON'

wof_df = pd.read_csv(os.path.join(input_data_dir,
                                  'Woff 2015.csv'),
                     encoding='latin-1')
wof_df.name = 'WOF'


# ### Create 'output' directory if it doesn't exist

# In[5]:


import os
if not os.path.exists(output_data_dir):
    os.makedirs(output_data_dir)


# ### Create columns for year, month, day, hour, and day type

# In[6]:


demand_nodes = [x for x in demand_df.columns if x != 'Datetime']

# Convert datetime to year, month, day, and hour
demand_df['Datetime'] = pd.to_datetime(demand_df['Datetime'])
demand_df['Year'] = demand_df['Datetime'].dt.strftime('%Y').astype(int)
demand_df['Month'] = demand_df['Datetime'].dt.strftime('%m').astype(int)
demand_df['Day'] = demand_df['Datetime'].dt.strftime('%d').astype(int)
demand_df['Hour'] = demand_df['Datetime'].dt.strftime('%H').astype(int)

# Create column for weekday/weekend
demand_df['Day-of-week'] = demand_df['Datetime'].dt.dayofweek
demand_df.loc[demand_df['Day-of-week'] < 5, 'Day-of-week'] = 'WD'
demand_df.loc[demand_df['Day-of-week'] != 'WD', 'Day-of-week'] = 'WE'


# ### Create dictionaries for 'seasons' and 'dayparts'

# In[7]:


seasons_dict = dict(zip(list(seasons_df['month']),
                        list(seasons_df['season'])
                       )
                   )

dayparts_dict = {i: [j, k] 
                 for i, j, k 
                 in zip(list(dayparts_df['daypart']), 
                        list(dayparts_df['start_hour']), 
                        list(dayparts_df['end_hour'])
                                             )
                }


# ### Create columns with 'seasons' and 'dayparts'

# In[8]:


demand_df['Season'] = demand_df['Month']
demand_df['Season'].replace(seasons_dict, inplace=True)

demand_df['Hour'] = demand_df['Hour'].map(lambda x: apply_timeshift(int(x), timeshift))
for daypart in dayparts_dict:
    if dayparts_dict[daypart][0] > dayparts_dict[daypart][1]: # loops over 24hrs
        demand_df.loc[(demand_df['Hour'] >= dayparts_dict[daypart][0]) |
                      (demand_df['Hour'] < dayparts_dict[daypart][1]),
                      'Daypart'] = daypart
    else:
        demand_df.loc[(demand_df['Hour'] >= dayparts_dict[daypart][0]) &
                  (demand_df['Hour'] < dayparts_dict[daypart][1]),
                  'Daypart'] = daypart


# ### Create column for timeslice with and without day-type

# In[9]:


if daytype_included:
    demand_df['TIMESLICE'] = (demand_df['Season'] +
                              demand_df['Day-of-week'] +
                              demand_df['Daypart'])
else:
    demand_df['TIMESLICE'] = (demand_df['Season'] +
                              demand_df['Daypart'])  


# ### Calculate YearSplit

# In[10]:


yearsplit = (demand_df['TIMESLICE'].
             value_counts(normalize = True).
             to_frame('VALUE').
             round(4).
             reset_index().
             rename({'index':'TIMESLICE'}, axis = 1))

yearsplit_final = pd.DataFrame(list(itertools.product(yearsplit['TIMESLICE'].unique(),
                                                      years)
                                   ),
                               columns = ['TIMESLICE', 'YEAR']
                              )
yearsplit_final = yearsplit_final.join(yearsplit.set_index('TIMESLICE'), 
                                       on = 'TIMESLICE')
yearsplit_final.to_csv(os.path.join(output_data_dir, 
                                    'YearSplit.csv'),
                       index=None)


# ### Calculate SpecifiedAnnualDemand and SpecifiedDemandProfile

# In[11]:


sp_demand_df = demand_df[[x 
                          for x in demand_df.columns 
                          if x in demand_nodes or
                          x == 'TIMESLICE']]
sp_demand_df = pd.melt(sp_demand_df,
                       id_vars = 'TIMESLICE',
                       value_vars = demand_nodes, 
                       var_name = 'node', 
                       value_name = 'demand')

sp_demand_df = sp_demand_df.groupby(['TIMESLICE', 'node'], 
                                    as_index = False).agg(sum)

# Calculate SpecifiedAnnualDemand
total_demand_df = sp_demand_df.groupby('node', 
                                       as_index = False).agg(sum)

total_demand_df.rename({'demand':'total_demand'}, 
                       axis = 1, 
                       inplace = True)

sp_demand_df = sp_demand_df.join(total_demand_df.set_index('node'), 
                                 on = 'node')

# Calculate SpecifiedDemandProfile

sp_demand_df['VALUE'] = (sp_demand_df['demand'] / 
                         sp_demand_df['total_demand'])

                
# Filter out country aggregate values for countries with multiple nodes 
country_with_nodes = list((sp_demand_df.loc[sp_demand_df['node'].str.len() > 6,
                                       'node'].
                      str[:-3].
                      unique()))

sp_demand_df = sp_demand_df.loc[~(sp_demand_df['node'].
                                  isin(country_with_nodes))]


# Rename COMMODITY based on naming convention. 
# Add 'XX' for countries without multiple nodes
sp_demand_df.loc[sp_demand_df['node'].str.len() <= 6, 
                 'FUEL'] = ('ELC' + 
                                 sp_demand_df['node'].
                                 str.split('-').
                                 str[1:].
                                 str.join("") +
                                 'XX02')

sp_demand_df.loc[sp_demand_df['node'].str.len() > 6, 
                 'FUEL'] = ('ELC' + 
                                 sp_demand_df['node'].
                                 str.split('-').
                                 str[1:].
                                 str.join("") +
                                 '02')

# Create master table for SpecifiedDemandProfile
sp_demand_df_final = pd.DataFrame(list(itertools.product(sp_demand_df['TIMESLICE'].unique(),
                                                         sp_demand_df['FUEL'].unique(),
                                                         years)
                                      ),
                                  columns = ['TIMESLICE', 'FUEL', 'YEAR']
                                 )
sp_demand_df_final = sp_demand_df_final.join(sp_demand_df.set_index(['TIMESLICE', 'FUEL']), 
                                             on = ['TIMESLICE', 'FUEL'])

# Add 'REGION' column and fill 'GLOBAL' throughout
sp_demand_df_final['REGION'] = 'GLOBAL'

total_demand_df_final = (sp_demand_df_final.
                         groupby(['REGION',
                                  'FUEL',
                                  'YEAR'],
                                 as_index = False)['total_demand'].
                         agg('mean').
                         rename({'total_demand':'VALUE'}, 
                                axis = 1)
                        )

# Convert SpecifiedAnnualDemand to required units
total_demand_df_final['VALUE'] = total_demand_df_final['VALUE'].mul(3.6*1e-6)

# Generate SpecifiedAnnualDemand.csv file 
#total_demand_df_final.to_csv(os.path.join(output_dir,'SpecifiedAnnualDemand.csv'), index=None)

# Generate SpecifiedDemandProfile.csv file 
sp_demand_df_final = sp_demand_df_final[['REGION',
                                         'FUEL',
                                         'TIMESLICE',
                                         'YEAR', 
                                         'VALUE']]

sp_demand_df_final.to_csv(os.path.join(output_data_dir,'SpecifiedDemandProfile.csv'), index=None)

# ### CapacityFactor

# In[12]:


datetime_ts_df = demand_df[['Datetime', 'TIMESLICE']]
capfac_all_df = pd.DataFrame(columns = ['REGION',
                                          'TECHNOLOGY',
                                          'TIMESLICE',
                                          'YEAR', 
                                          'VALUE'])

def capacity_factor(df):
    df['Datetime'] = pd.to_datetime(df['Datetime'])
    capfac_df = (df.
                 set_index('Datetime').
                 join(datetime_ts_df.
                      set_index('Datetime'),
                      on = 'Datetime')
                )
    capfac_nodes = [x for 
                    x in 
                    capfac_df.columns 
                    if x 
                    not in ['Datetime', 'TIMESLICE']]
    capfac_df = capfac_df.reset_index().drop('Datetime', 
                                             axis = 1)
    capfac_df = pd.melt(capfac_df,
                        id_vars = 'TIMESLICE',
                        value_vars = capfac_nodes, 
                        var_name = 'node', 
                        value_name = 'VALUE')
    capfac_df = (capfac_df.
                 groupby(['TIMESLICE', 'node'],
                         as_index = False).
                 agg('mean')
                )
    capfac_df['VALUE'] = (capfac_df['VALUE'].
                          div(100).
                          round(4)
                         )
    
    ## Filter out country aggregate values for countries with multiple nodes 
    capfac_df = capfac_df.loc[~(capfac_df['node'].
                                isin(country_with_nodes))]

    # Rename COMMODITY based on naming convention. 
    # Add 'XX' for countries without multiple nodes
    capfac_df.loc[capfac_df['node'].str.len() <= 6, 
                  'TECHNOLOGY'] = ('PWR' +
                                   df.name +
                                   capfac_df['node'].
                                   str.split('-').
                                   str[1:].
                                   str.join("") +
                                   'XX01')
    
    capfac_df.loc[capfac_df['node'].str.len() > 6, 
                  'TECHNOLOGY'] = ('PWR' + 
                                   df.name +
                                   capfac_df['node'].
                                   str.split('-').
                                   str[1:].
                                   str.join("") +
                                   '01')
    
    # Create master table for CapacityFactor
    capfac_df_final = pd.DataFrame(list(itertools.product(capfac_df['TIMESLICE'].unique(),
                                                          capfac_df['TECHNOLOGY'].unique(),
                                                          years)
                                      ),
                                   columns = ['TIMESLICE', 'TECHNOLOGY', 'YEAR']
                                 )
    capfac_df_final = (capfac_df_final.
                       join(capfac_df.
                            set_index(['TIMESLICE', 'TECHNOLOGY']), 
                            on = ['TIMESLICE', 'TECHNOLOGY'])
                      )
    
    # Add 'REGION' column and fill 'GLOBAL' throughout
    capfac_df_final['REGION'] = 'GLOBAL'
    
    capfac_df_final = capfac_df_final[['REGION',
                                       'TECHNOLOGY',
                                       'TIMESLICE',
                                       'YEAR', 
                                       'VALUE']]
    
    return capfac_df_final


for each in [hyd_df_processed, csp_df, spv_df, won_df, wof_df]:
    capfac_all_df = capfac_all_df.append(capacity_factor(each),
                                         ignore_index = True)
    
capfac_all_df.to_csv(os.path.join(output_data_dir, 
                                  'CapacityFactor.csv'),
                     index=None)


# ## Create csv for TIMESLICE 

# In[13]:


time_slice_list = list(demand_df['TIMESLICE'].unique())
time_slice_df = pd.DataFrame(time_slice_list, columns = ['VALUE'])
time_slice_df.to_csv(os.path.join(output_data_dir, 
                                  'TIMESLICE.csv'),
                     index=None)

logging.info('Time Slicing Completed')
