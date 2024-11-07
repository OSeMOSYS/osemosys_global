"""Functions to set ReserveMarginTagTechnology."""
import pandas as pd
import itertools

from data import get_years

def set_reserve_margin_technologies(margins_technologies, 
                                    tech_set, start_year, 
                                    end_year, region_name):
    
    years = get_years(start_year, end_year)
    data = pd.DataFrame()
    
    '''Check if transmission is added to the config to contribute to
    reserves and develop outputs.'''
    for tech, val in margins_technologies.items():
        if tech == 'TRN':
            rm_techs = [
                x
                for x in tech_set["VALUE"].unique()
                if x.startswith("TRN")
                ]
    
    '''Check if generator/storage tech is added to the config to 
    contribute to reserves and develop outputs.'''
        else:
            rm_techs = [
                x
                for x in tech_set["VALUE"].unique()
                if x.startswith("PWR")
                if x[3:6] == tech
            ]
        
        tech_data = pd.DataFrame(
            list(itertools.product([region_name], rm_techs, years)),
            columns=["REGION", "TECHNOLOGY", "YEAR"],
        )
        
        tech_data['VALUE'] = val / 100
        data = pd.concat([data, tech_data])

    return data