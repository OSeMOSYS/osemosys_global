"""Functions to set operational life."""

import pandas as pd

from utils import apply_dtypes

def set_op_life_storage(storage_set, op_life_dict, 
                        op_life_base, region_name):

    # Create Operational Life data
    op_life_techs = list(storage_set['VALUE'].unique())
    
    op_life_storage_out = []
    op_life_out = []
    
    # transmission technologies 
    for op_life_tech in op_life_techs:
        op_life_storage_out.append([region_name, op_life_tech, op_life_dict.get(op_life_tech[:3])])
        op_life_out.append([region_name, 'PWR' + op_life_tech, op_life_dict.get(op_life_tech[:3])])
        
    op_life_storage_final = pd.DataFrame(op_life_storage_out, columns = ['REGION', 'STORAGE', 'VALUE'])
    
    op_life_final = pd.DataFrame(op_life_out, columns = ['REGION', 'TECHNOLOGY', 'VALUE'])
    op_life_final = pd.concat([op_life_base, op_life_final])
    
    op_life_storage_final = apply_dtypes(op_life_storage_final, "OperationalLife")
    op_life_final = apply_dtypes(op_life_final, "OperationalLife")
    
    return op_life_final, op_life_storage_final