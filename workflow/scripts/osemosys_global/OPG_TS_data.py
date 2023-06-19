#!/usr/bin/env python
# coding: utf-8

# # OSeMOSYS-PLEXOS global model: TS-dependent parameters

# ### Import modules


import pandas as pd
import datetime
import numpy as np
import itertools
import seaborn as sns; sns.set()
import matplotlib
import matplotlib.pyplot as plt
import urllib
import os
import time
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
custom_nodes_dir = config_paths.custom_nodes_dir
geographic_scope = config.get('geographic_scope')
seasons = config.get('seasons')
daytype = config.get('daytype')
dayparts = config.get('dayparts')
reserve_margin = config.get('reserve_margin')

# Check for custom nodes directory
try:
    os.makedirs(custom_nodes_dir)
except FileExistsError:
    pass

region_name = config.region_name
custom_nodes = config.get('nodes_to_add')

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

# Read renewable profile files
csp_df = pd.read_csv(os.path.join(input_data_dir,
                                  'CSP 2015.csv'),
                     encoding='latin-1')
csp_df.name = 'CSP'

spv_df = pd.read_csv(os.path.join(input_data_dir,
                                  'SolarPV 2015.csv'),
                     encoding='latin-1')
if custom_nodes:
    spv_df_custom = pd.read_csv(os.path.join(custom_nodes_dir,
                                             'RE_profiles_SPV.csv'),
                                             encoding='latin-1')
    spv_df_custom.drop(['Datetime'],
                        axis=1,
                        inplace=True)
    spv_df = pd.concat([spv_df, spv_df_custom], axis=1)

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
if custom_nodes:
    hyd_df_custom = pd.read_csv(os.path.join(custom_nodes_dir,
                                             'RE_profiles_HYD.csv'),
                                             encoding='latin-1')
    hyd_df = pd.concat([hyd_df, hyd_df_custom])
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
if custom_nodes:
    won_df_custom = pd.read_csv(os.path.join(custom_nodes_dir,
                                             'RE_profiles_WON.csv'),
                                             encoding='latin-1')
    won_df_custom.drop(['Datetime'],
                        axis=1,
                        inplace=True)
    won_df = pd.concat([won_df, won_df_custom], axis=1)
won_df.name = 'WON'

wof_df = pd.read_csv(os.path.join(input_data_dir,
                                  'Woff 2015.csv'),
                     encoding='latin-1')
if custom_nodes:
    wof_df_custom = pd.read_csv(os.path.join(custom_nodes_dir,
                                             'RE_profiles_WOF.csv'),
                                             encoding='latin-1')
    wof_df_custom.drop(['Datetime'],
                        axis=1,
                        inplace=True)
    wof_df = pd.concat([wof_df, wof_df_custom], axis=1)
wof_df.name = 'WOF'


# ### Create 'output' directory if it doesn't exist

import os
if not os.path.exists(output_data_dir):
    os.makedirs(output_data_dir)


# ### Create columns for year, month, day, hour, and day type

if custom_nodes:
    demand_nodes = [x for x in demand_df.columns if x != 'Datetime'] + custom_nodes
else:
    demand_nodes = [x for x in demand_df.columns if x != 'Datetime']
# Convert datetime to year, month, day, and hour
demand_df['Datetime'] = pd.to_datetime(demand_df['Datetime'])
demand_df['Year'] = demand_df['Datetime'].dt.strftime('%Y').astype(int)
demand_df['Month'] = demand_df['Datetime'].dt.strftime('%m').astype(int)
demand_df['Day'] = demand_df['Datetime'].dt.strftime('%d').astype(int)
demand_df['Hour'] = demand_df['Datetime'].dt.strftime('%H').astype(int)

if custom_nodes:
    custom_sp_demand_profile = pd.read_csv(os.path.join(input_data_dir,
                                           "custom_nodes",
                                           "specified_demand_profile.csv"))
    demand_df = pd.merge(demand_df,
                         custom_sp_demand_profile,
                         how='left',
                         on=['Month','Day','Hour'])

# Create column for weekday/weekend
demand_df['Day-of-week'] = demand_df['Datetime'].dt.dayofweek
demand_df.loc[demand_df['Day-of-week'] < 5, 'Day-of-week'] = 'WD'
demand_df.loc[demand_df['Day-of-week'] != 'WD', 'Day-of-week'] = 'WE'


# ### Create dictionaries for 'seasons' and 'dayparts'

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

if daytype_included:
    demand_df['TIMESLICE'] = (demand_df['Season'] +
                              demand_df['Day-of-week'] +
                              demand_df['Daypart'])
else:
    demand_df['TIMESLICE'] = (demand_df['Season'] +
                              demand_df['Daypart'])  


# ### Calculate YearSplit

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
sp_demand_df.loc[sp_demand_df['node'].str.len() == 5, 
                 'FUEL'] = ('ELC' + 
                            sp_demand_df['node'] +
                            '02')

sp_demand_df.loc[sp_demand_df['node'].str.len() == 6, 
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
sp_demand_df_final['VALUE'] = sp_demand_df_final['VALUE'].round(2)
sp_demand_df_final = sp_demand_df_final[['REGION',
                                         'FUEL',
                                         'TIMESLICE',
                                         'YEAR', 
                                         'VALUE']]
sp_demand_df_final.drop_duplicates(subset=['REGION','TIMESLICE','FUEL','YEAR'],
                                   keep='last',
                                   inplace=True)
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

capfac_all_df.drop_duplicates(subset=['REGION','TECHNOLOGY','TIMESLICE','YEAR'],
                              keep='last',
                              inplace=True)    
capfac_all_df.to_csv(os.path.join(output_data_dir, 
                                  'CapacityFactor.csv'),
                     index=None)

# ## Create csv for TIMESLICE
time_slice_list = list(demand_df['TIMESLICE'].unique())
time_slice_df = pd.DataFrame(time_slice_list, columns = ['VALUE'])
time_slice_df.to_csv(os.path.join(output_data_dir, 
                                  'TIMESLICE.csv'),
                     index=None)

'''
def add_storage(region_name, 
                years, 
                output_data_dir, 
                demand_nodes, 
                time_slice_list,
                seasons,
                daytype,
                dayparts):
                
            '''

demand_nodes = list(set(list(sp_demand_df_final['FUEL'].str[3:8])))

# Create SET STORAGE
storage_set = [('BAT' + x +'01') for x in demand_nodes 
                if x[:3] in geographic_scope]
df_storage_set = pd.DataFrame(storage_set,
                              columns=['VALUE'])
df_storage_set.to_csv(os.path.join(output_data_dir,
                                    'STORAGE.csv'),
                        index=None)
# Add storage technologies to SET TECHNOLOGY
storage_techs = [('PWRBAT' + x +'01') for x in demand_nodes 
                if x[:3] in geographic_scope]
df_storage_techs = pd.DataFrame(storage_techs,
                                columns=['VALUE'])

wait_time = 0
while not os.path.exists(os.path.join(output_data_dir, 'TECHNOLOGY.csv')):
    time.sleep(5)
    wait_time += 1
    if wait_time > 20 : break

set_techonology = pd.read_csv(os.path.join(output_data_dir,
                                           'TECHNOLOGY.csv'))
set_technology = pd.concat([set_techonology, df_storage_techs])
set_technology.to_csv(os.path.join(output_data_dir,
                                    'TECHNOLOGY.csv'),
                        index=None)

# Add InputActivityRatio and OutputActivityRatio
# InputActivityRatio
df_storage_iar = pd.DataFrame(list(itertools.product([region_name],
                                                     storage_techs,
                                                     years,
                                                     [1])),
                              columns=['REGION',
                                       'TECHNOLOGY',
                                       'YEAR',
                                       'MODE_OF_OPERATION']
                              )
df_storage_iar['VALUE'] = 1
df_storage_iar['FUEL'] = 'ELC' + df_storage_iar['TECHNOLOGY'].str[6:11] + '01'
df_storage_iar = df_storage_iar[['REGION',
                                 'TECHNOLOGY',
                                 'FUEL',
                                 'MODE_OF_OPERATION',
                                 'YEAR',
                                 'VALUE']]

wait_time = 0
while not os.path.exists(os.path.join(output_data_dir, 'InputActivityRatio.csv')):
    time.sleep(5)
    wait_time += 1
    if wait_time > 20 : break
df_iar = pd.read_csv(os.path.join(output_data_dir,
                                  'InputActivityRatio.csv'))
df_iar = pd.concat([df_iar, df_storage_iar])
df_iar.to_csv(os.path.join(output_data_dir,
                           'InputActivityRatio.csv'),
              index=None)
time.sleep(20)

# OutputActivityRatio
df_storage_oar = pd.DataFrame(list(itertools.product([region_name],
                                                     storage_techs,
                                                     years,
                                                     [2])),
                              columns=['REGION',
                                       'TECHNOLOGY',
                                       'YEAR',
                                       'MODE_OF_OPERATION']
                              )
df_storage_oar['VALUE'] = 1
df_storage_oar['FUEL'] = 'ELC' + df_storage_oar['TECHNOLOGY'].str[6:11] + '01'
df_storage_oar = df_storage_oar[['REGION',
                                 'TECHNOLOGY',
                                 'FUEL',
                                 'MODE_OF_OPERATION',
                                 'YEAR',
                                 'VALUE']]

wait_time = 0
while not os.path.exists(os.path.join(output_data_dir, 'OutputActivityRatio.csv')):
    time.sleep(5)
    wait_time += 1
    if wait_time > 20 : break
df_oar = pd.read_csv(os.path.join(output_data_dir,
                                  'OutputActivityRatio.csv'))
df_oar = pd.concat([df_oar, df_storage_oar])
df_oar.to_csv(os.path.join(output_data_dir,
                           'OutputActivityRatio.csv'),
              index=None)
time.sleep(20)

# Create TechnologyToStorage and TechnologyFromStorage

df_tech_storage = pd.DataFrame(columns=['REGION',
                                        'TECHNOLOGY',
                                        'STORAGE',
                                        'MODE_OF_OPERATION'])

for each_node in [x for x in demand_nodes if x[:3] in geographic_scope]:
    df_ts_temp = pd.DataFrame(list(itertools.product([region_name],
                                                        ['PWRBAT' + each_node +'01'],
                                                        ['BAT' + each_node +'01'],
                                                        [1,2])
                                    ),
                                columns=['REGION',
                                        'TECHNOLOGY',
                                        'STORAGE',
                                        'MODE_OF_OPERATION']
                                )
    df_tech_storage = pd.concat([df_tech_storage, df_ts_temp])

df_ttos = df_tech_storage.copy()
df_tfroms = df_tech_storage.copy()


# TechnologyToStorage

df_ttos.loc[df_ttos['MODE_OF_OPERATION'] == 1,
            'VALUE'] = 1.0
df_ttos.loc[df_ttos['MODE_OF_OPERATION'] == 2,
            'VALUE'] = 0.0
df_ttos['VALUE'] = df_ttos['VALUE'].astype(float)
df_ttos.to_csv(os.path.join(output_data_dir,
                            'TechnologyToStorage.csv'),
               index=None)

# TechnologyFromStorage

df_tfroms.loc[df_tfroms['MODE_OF_OPERATION'] == 1,
                'VALUE'] = 0.0
df_tfroms.loc[df_tfroms['MODE_OF_OPERATION'] == 2,
                'VALUE'] = 1.0
df_tfroms['VALUE'] = df_tfroms['VALUE'].astype(float)
df_tfroms.to_csv(os.path.join(output_data_dir,
                                'TechnologyFromStorage.csv'),
                 index=None)
    
# Create Conversionls, Conversionld, and Conversionlh

# Conversionls
df_ls = pd.DataFrame(list(itertools.product(time_slice_list,
                                            list(range(1,len(seasons)+1))
                                            )
                            ),
                        columns=['TIMESLICE',
                                'SEASON']
                        )
df_ls.loc[df_ls['TIMESLICE'].str[1:2].astype(int) == df_ls['SEASON'],
            'VALUE'] = 1
df_ls.fillna(0, inplace=True)
df_ls.to_csv(os.path.join(output_data_dir,
                            'Conversionls.csv'),
             index=None)

df_season_set = pd.DataFrame(list(range(1,len(seasons)+1)),
                             columns=['VALUE'])
df_season_set.to_csv(os.path.join(output_data_dir,
                                  'SEASON.csv'),
                       index=None)

# Conversionld
df_ld = pd.DataFrame(list(itertools.product(time_slice_list,
                                            [1]
                                            )
                            ),
                        columns=['TIMESLICE',
                                'DAYTYPE']
                        )
df_ld['VALUE'] = 1
df_ld.fillna(0, inplace=True)
df_ld.to_csv(os.path.join(output_data_dir,
                            'Conversionld.csv'),
             index=None)
df_daytype_set = pd.DataFrame([1],
                              columns=['VALUE'])
df_daytype_set.to_csv(os.path.join(output_data_dir,
                                   'DAYTYPE.csv'),
                       index=None)

# Conversionlh
df_lh = pd.DataFrame(list(itertools.product(time_slice_list,
                                            list(range(1,len(dayparts)+1))
                                                    )
                            ),
                        columns=['TIMESLICE',
                                'DAILYTIMEBRACKET']
                        )
df_lh.loc[df_lh['TIMESLICE'].str[3:].astype(int) == df_lh['DAILYTIMEBRACKET'],
            'VALUE'] = 1
df_lh.fillna(0, inplace=True)
df_lh.to_csv(os.path.join(output_data_dir,
                            'Conversionlh.csv'),
             index=None)
df_dayparts_set = pd.DataFrame(list(range(1,len(dayparts)+1)),
                               columns=['VALUE'])
df_dayparts_set.to_csv(os.path.join(output_data_dir,
                                    'DAILYTIMEBRACKET.csv'),
                       index=None)

# Daysplit

daysplit = {}
for dp, hr in dayparts_raw.items():
    daysplit[int(dp[1:])] = (hr[1] - hr[0])/8760

df_daysplit = pd.DataFrame(itertools.product(list(range(1,len(dayparts)+1)),
                                             years),
                           columns=['DAILYTIMEBRACKET',
                                    'YEAR'])
df_daysplit['VALUE'] = df_daysplit['DAILYTIMEBRACKET'].map(daysplit)
df_daysplit = df_daysplit[['DAILYTIMEBRACKET',
                           'YEAR',
                           'VALUE']]
df_daysplit['VALUE'] = df_daysplit['VALUE'].round(4)
df_daysplit.to_csv(os.path.join(output_data_dir,
                                'DaySplit.csv'),
                   index=None)

# CapitalCostStorage
storage_set = [('BAT' + x +'01') for x in demand_nodes 
                if x[:3] in geographic_scope]
df_cap_cost_storage = pd.DataFrame(list(itertools.product(storage_set,
                                                          years)
                                        ),
                                   columns=['STORAGE',
                                            'YEAR']
                                   )
df_cap_cost_storage['STORAGE_TYPE'] = df_cap_cost_storage['STORAGE'].str[:3]
storage_costs = pd.read_csv(os.path.join(input_data_dir,
                                         'storage_costs.csv'))

storage_costs_df = pd.DataFrame(list(itertools.product(storage_costs['STORAGE_TYPE'].unique(),
                                                       list(range(storage_costs['YEAR'].min(),
                                                                  storage_costs['YEAR'].max()+1)))
                                        ),
                                   columns=['STORAGE_TYPE',
                                            'YEAR']
                                   )
storage_costs_df = storage_costs_df.merge(storage_costs,
                                          how='left',
                                          on=['STORAGE_TYPE', 'YEAR'])
storage_costs_df = storage_costs_df.interpolate()
df_cap_cost_storage = df_cap_cost_storage.merge(storage_costs_df,
                                                how='left',
                                                on=['STORAGE_TYPE', 'YEAR'])
df_cap_cost_storage['VALUE'] = df_cap_cost_storage['VALUE'].mul(1e6/3600)
df_cap_cost_storage['REGION'] = region_name
df_cap_cost_storage = df_cap_cost_storage[['REGION',
                                           'STORAGE',
                                           'YEAR',
                                           'VALUE']]
df_cap_cost_storage.to_csv(os.path.join(output_data_dir,
                                        'CapitalCostStorage.csv'),
                           index=None)

# CapacityToActivityUnit for Storage


# ReserveMargin
    
df_rm = pd.DataFrame(years,
                     columns=['YEAR'])
for rm, rm_params in reserve_margin.items():
    df_rm.loc[df_rm['YEAR'].between(rm_params[1], rm_params[2]),
                'VALUE'] = (1 + rm_params[0]/100)

df_rm = df_rm.interpolate()    
df_rm['REGION'] = region_name
df_rm = df_rm[['REGION',
               'YEAR',
               'VALUE']]
df_rm.to_csv(os.path.join(output_data_dir,
                            'ReserveMargin.csv'),
             index=None)

# ReserveMarginTagTechnology
df_rmtt = pd.read_csv(os.path.join(output_data_dir,
                                    'TECHNOLOGY.csv'))
reserve_margin_techs = ['COA',
                        'COG',
                        'OCG',
                        'CCG',
                        'PET',
                        'URN',
                        'OIL',
                        'OTH',
                        'BIO',
                        'HYD',
                        'GEO',
                        'SPV',
                        'WON'
                        ]
rm_techs = [x for x in df_rmtt['VALUE'].unique()
            if x.startswith('PWR')
            if x[3:6] in reserve_margin_techs]      
df_rmtt = pd.DataFrame(list(itertools.product([region_name],
                                              rm_techs,
                                              years,
                                              [1])),
                       columns=['REGION',
                                'TECHNOLOGY',
                                'YEAR',
                                'VALUE']
                       )                  
df_rmtt.to_csv(os.path.join(output_data_dir,
                            'ReserveMarginTagTechnology.csv'),
               index=None)

# ReserveMarginTagFuel
df_rmtf = pd.read_csv(os.path.join(output_data_dir,
                                   'FUEL.csv'))
rm_fuels = [x for x in df_rmtf['VALUE'].unique()
            if x.startswith('ELC')
            if x.endswith('01')]
df_rmtf = pd.DataFrame(list(itertools.product([region_name],
                                              rm_fuels,
                                              years,
                                              [1])),
                       columns=['REGION',
                                'FUEL',
                                'YEAR',
                                'VALUE']
                       )    
df_rmtf.to_csv(os.path.join(output_data_dir,
                            'ReserveMarginTagFuel.csv'),
               index=None)
logging.info('Time Slicing Completed')
