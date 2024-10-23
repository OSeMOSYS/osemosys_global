"""Functions to set operational life."""

import pandas as pd

from utils import apply_dtypes

def set_op_life_transmission(df_oar_final, op_life_dict, 
                             op_life_base, region_name):

    # Create Operational Life data
    op_life_techs = list(df_oar_final['TECHNOLOGY'].unique())
    
    op_life_trn = [tech for tech in op_life_techs if tech[0:3] == 'TRN']
    
    op_life_out = []
    
    # transmission technologies 
    for op_life_tech in op_life_trn:
        op_life_out.append([region_name, op_life_tech, op_life_dict.get('TRN')])
        
    op_life_trn_final = pd.DataFrame(op_life_out, columns = ['REGION', 
                                                             'TECHNOLOGY', 'VALUE'])

    op_life_trn_final = apply_dtypes(op_life_trn_final, "OperationalLife")
    
    op_life_trn_final = pd.concat([op_life_base, op_life_trn_final])

    return op_life_trn_final