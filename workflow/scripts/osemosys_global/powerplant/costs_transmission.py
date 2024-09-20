"""Function to calculate residual capacity for powerplant technologies."""

from data_transmission import format_transmission_name

from constants import(
    region_name,
    years
    )

def get_transmission_costs(df_trn_lines, df_oar_final):
    '''Gets electrical transmission capital and fixed cost per technology. 

    Both the capital costs and fixed cost are written out to avoid having 
    to read in the excel file twice 
    
    Returns: 
        df_trans_capex: DataFrame with the columns 'TECHNOLOGY' and 'VALUE' 
            representing CAPITAL cost in millions of dollars per year. 
        df_trans_fix: DataFrame with the columns 'TECHNOLOGY' and 'VALUE' 
            representing FIXED cost in millions of dollars per year. 
    '''

    # Drop unneeded columns
    df = df_trn_lines.drop(
        [
            "Line",
            "KM distance",
            "HVAC/HVDC/Subsea",
            "Losses",
            "Unnamed: 8",
            "Line Max Size (MW)",
            "Unnamed: 10",
            "Unnamed: 11",
            "Unnamed: 12",
            "Subsea lines",
            "Unnamed: 14"
        ],
        axis=1,
    )

    # Use to/from codes to create a TECHNOLOGY columns
    df = format_transmission_name(df)

    # Changes units
    # Raw given in dollars -> model units in million dollars
    df = df.rename(columns={'Annual FO&M (3.5% of CAPEX) ($2010 in $000)':'O&M'})
    df = df.rename(columns={'Build Cost ($2010 in $000)':'CAPEX'})
    df['O&M'] = df['O&M'].astype(float) / 1000
    df['CAPEX'] = df['CAPEX'].astype(float) / 1000
    df = df.round({'O&M': 3, 'CAPEX': 3, })

    # Separate out fixed and capex and return values 
    df_trans_capex = df.drop(['O&M'], axis=1)
    df_trans_capex = df_trans_capex.rename(columns={'CAPEX':'VALUE'})
    df_trans_fix = df.drop(['CAPEX'], axis=1)
    df_trans_fix = df_trans_fix.rename(columns={'O&M':'VALUE'})

    df_trans_capex['REGION'] = region_name
    df_trans_capex['YEAR'] = [years] * len(df_trans_capex)
    df_trans_capex = df_trans_capex.explode('YEAR')

    df_trans_fix['REGION'] = region_name
    df_trans_fix['YEAR'] = [years] * len(df_trans_fix)
    df_trans_fix = df_trans_fix.explode('YEAR')

    # Filter out techs that don't have activity ratios 
    df_trans_capex = df_trans_capex.loc[
        df_trans_capex['TECHNOLOGY'].isin(df_oar_final['TECHNOLOGY'])]
    df_trans_fix = df_trans_fix.loc[
        df_trans_fix['TECHNOLOGY'].isin(df_oar_final['TECHNOLOGY'])]
    
    return df_trans_capex, df_trans_fix