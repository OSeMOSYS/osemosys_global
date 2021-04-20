#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import os
import shutil


# In[ ]:

input_dir = r'../../data/'
output_dir = r'../../osemosys_global_model/data/'                       

def create_csv_files(path):
    csv_file_dict = {} 
    for each_csv in os.listdir(r'../../simplicity/data/'):  
        if each_csv not in os.listdir(path):
            #copy default values csv with data
            if each_csv == "default_values.csv":
                shutil.copy(r'../../simplicity/data/default_values.csv',path)
            #copy file names and column headers only
            else:
                csv_df_in = pd.read_csv(os.path.join(r'../../simplicity/data/',
                                                  each_csv))
                csv_df_out = pd.DataFrame(columns = list(csv_df_in.columns))
                csv_df_out.to_csv(os.path.join(path, 
                                               each_csv),
                                 index = None)
                csv_file_dict[each_csv] = list(csv_df_in.columns)
    
    return csv_file_dict


# In[ ]:


create_csv_files(output_dir)


# In[ ]:




