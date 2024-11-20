"""Function to set AnnualEmissionLimit."""

import pandas as pd
import itertools

from data import get_years

def add_emission_limits(emissions, emission_limit, ember, 
                        start_year, end_year, region_name):
    
    # Set empty AnnualEmissionLimit df.
    annual_emission_limit = pd.DataFrame(
        columns=["REGION", "EMISSION", "YEAR", "VALUE"]
    ) 
    
    # Check if emission limits exist.
    if emission_limit:
        emissions_list = list(emissions["VALUE"])
        years = get_years(start_year, end_year)
        
        # Order emission limits based on emission type, country and year entries.
        emission_limit_ordered = pd.DataFrame()
        
        for el_params in emission_limit:
            param = pd.DataFrame([el_params], columns = ['EMISSION', 'COUNTRY', 'TYPE', 
                                                         'YEAR', 'VALUE'])
            if emission_limit_ordered.empty:
                emission_limit_ordered = param
                
            else:
                emission_limit_ordered = pd.concat([emission_limit_ordered, param])
            
        emission_limit_ordered = emission_limit_ordered.sort_values(
            by = ['EMISSION', 'COUNTRY', 'YEAR']).reset_index(drop = True)

        # Create dataframe template for adding limits.
        template = pd.DataFrame(
            list(itertools.product(emissions_list, years)), columns=["EMISSION", "YEAR"]
        )

        # Loop through emissions       
        for emission in emission_limit_ordered['EMISSION'].unique():
            # Loop through countries
            for country in emission_limit_ordered['COUNTRY'].unique():
                prev_data_year = None
                prev_data_value = None

                # Set input limits for given country and emission
                limits = emission_limit_ordered.loc[(emission_limit_ordered[
                    'COUNTRY'] == country) & (emission_limit_ordered[
                        'EMISSION'] == emission)]
                
                # Filter template df for given country and emission
                data = template.copy().loc[template["EMISSION"].isin(
                    [emission + country])]
                     
                # Loop through limits per country
                for idx, row in limits.iterrows():
                    data_type = row.TYPE
                    data_year = row.YEAR
                    data_value = row.VALUE
                    
                    # Check if a limit for earlier years is set
                    if not prev_data_year:
                        # If not check if limit type is LINEAR. If TRUE, set historic
                        # emission values to be able to set a linear relationship in absence
                        # of previous year entries.
                        if data_type == 'LINEAR':
                            ember_data = ember.copy().loc[
                                (ember["EMISSION"].isin([emission + country]))]
                            # Check if year for country specific EMBER data exists in horizon.
                            for year in ember_data["YEAR"].unique():
                                if year in data["YEAR"].values:
                                    # Set baseline year data for LINEAR
                                    data.loc[(data["YEAR"] == year),"VALUE"] = ember_data.loc[
                                        (ember_data["YEAR"] == year)]['VALUE'].iloc[0]

                        # Set baseline year data for POINT
                        data.loc[(data["YEAR"] == data_year),"VALUE"] = data_value     
                    
                    else:
                        # If previous targets exist check if limit type is LINEAR. 
                        # If TRUE, only set target for target data year. If FALSE,
                        # set VALUE for all years between previous target year + 1 and 
                        # the current target year.
                        if data_type == 'LINEAR':
                            data.loc[(data["YEAR"] == data_year),"VALUE"] = data_value
                        else:
                            data.loc[data["YEAR"].between(prev_data_year, data_year)
                                     ,"VALUE"] = prev_data_value
                            
                            data.loc[(data["YEAR"] == data_year),"VALUE"] = data_value
                  
                    prev_data_year = data_year
                    prev_data_value = data_value
                    
                # Once done with a given country reset previous year limits.
                prev_data_year = None
                prev_data_value = None
                
                # Interpolate values for LINEAR targets.
                data['VALUE'] = data['VALUE'].interpolate()
      
                # Format df
                data.dropna(axis=0, inplace=True)
                data["REGION"] = region_name
                data = data[["REGION", "EMISSION", "YEAR", "VALUE"]]
                
                if annual_emission_limit.empty:
                    annual_emission_limit = data
                    
                else:
                    annual_emission_limit = pd.concat([annual_emission_limit, data])

    return annual_emission_limit