"""Function to calculate residual capacity for transmission technologies."""
import pandas as pd

from data import get_years

def res_capacity_transmission(df_exist_corrected, df_plan_corrected, 
                              transmission_existing, transmission_planned,
                              res_cap_base, op_life_dict, 
                              start_year, end_year, region_name,
                              retirement_year_transmission, 
                              planned_build_year_transmission):
    
    df_res_cap = None
    
    if transmission_existing:
        """Set build year for existing capacities as the differential between assumed 
        retirement year and the set operational life for transmission technologies."""
        df_res_cap_exist = df_exist_corrected.copy()
        df_res_cap_exist['YEAR'] = retirement_year_transmission - op_life_dict.get('TRN')
        
    if transmission_planned:
        """Set assumed build year for planned capacities from the GTD dataset that do not have 
        a commissioning year attached."""
        df_res_cap_plan = df_plan_corrected.copy()
        df_res_cap_plan.loc[df_res_cap_plan['YEAR'] == '-', 'YEAR'] = planned_build_year_transmission
    
    if transmission_existing and transmission_planned:
        # Combine dfs and filter out rows without capacities.
        df_res_cap = pd.concat([df_res_cap_exist, df_res_cap_plan])
    if transmission_existing and not transmission_planned:
        df_res_cap = df_res_cap_exist.copy()
    if not transmission_existing and transmission_planned:
        df_res_cap = df_res_cap_plan.copy()
              
    if df_res_cap is not None:
        df_res_cap = df_res_cap.loc[df_res_cap['VALUE'] != 0
                                              ].reset_index(drop = True
                                                            ).rename(columns = {'YEAR' :'build_year'})
                                                            
        # Set retirement years.                                                    
        df_res_cap['build_year'] = df_res_cap['build_year'].astype(int)
        df_res_cap['retirement_year'] = df_res_cap['build_year'] + op_life_dict.get('TRN')
        
        # Set residual capacity for all model horizon years.
        df_res_cap['YEAR'] = [get_years(start_year, end_year)] * len(df_res_cap)
        df_res_cap = df_res_cap.explode('YEAR').reset_index(drop = True)
        
        # Convert from MW to GW.
        df_res_cap['VALUE'] = df_res_cap['VALUE'] / 1000
        
        df_res_cap.loc[~
            ((df_res_cap["YEAR"] >= df_res_cap["build_year"])
            & (df_res_cap["YEAR"] <= df_res_cap["retirement_year"])),
            "VALUE",
        ] = 0
        
        # Group capacities by year and technology.
        df_res_cap = df_res_cap.groupby(['TECHNOLOGY', 'YEAR']).sum('VALUE').reset_index()
        
        # Reorder columns
        df_res_cap['REGION'] = region_name
        df_res_cap = df_res_cap[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]
        
        # Combine powerplant and transmission resisdual capacity df's.
        df_res_cap = pd.concat([res_cap_base, df_res_cap]).reset_index(drop = True)
        
    else:
        df_res_cap = res_cap_base.copy()

    return df_res_cap