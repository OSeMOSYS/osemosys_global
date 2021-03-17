#!/usr/bin/env python
# coding: utf-8

# # Filter osemosys_global datapackaged based on user-defined geographic scope

# ### Import modules

# In[ ]:


import pandas as pd
import os


# In[ ]:

input_dir = r'../../data/'
output_dir = r'../../osemosys_global_model/data/'                     

geographic_scope = ['IND', 'CHN']
model_name = 'geo_filter_test_v4'

if not os.path.exists(os.path.join(os.getcwd(),
                           model_name,
                           'data')):
    os.makedirs(os.path.join(os.getcwd(),
                           model_name,
                           'data'))

for each_csv in (os.listdir(os.path.join(os.getcwd(),
                                         output_dir))):
    
    df = pd.read_csv(os.path.join(os.getcwd(),
                                  output_dir,
                                  each_csv)
                    )
    if not df.empty:
        if 'TECHNOLOGY' in df.columns:
            df = df.loc[df['TECHNOLOGY'].str[3:6].isin(geographic_scope) | 
                        df['TECHNOLOGY'].str[6:9].isin(geographic_scope) | 
                        df['TECHNOLOGY'].str[8:11].isin(geographic_scope)]
            
        if 'FUEL' in df.columns:
            df = df.loc[df['FUEL'].str[3:6].isin(geographic_scope) | 
                        df['FUEL'].str[6:9].isin(geographic_scope)]
        
        if each_csv == 'FUEL.csv':
            df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                        df['VALUE'].str[6:9].isin(geographic_scope)]
        
        if each_csv == 'TECHNOLOGY.csv':
            df = df.loc[df['VALUE'].str[3:6].isin(geographic_scope) | 
                        df['VALUE'].str[6:9].isin(geographic_scope) | 
                        df['VALUE'].str[8:11].isin(geographic_scope)]
        
    
    df.to_csv(os.path.join(os.getcwd(),
                           model_name,
                           'data',
                           each_csv),
              index = None)


# In[ ]:




