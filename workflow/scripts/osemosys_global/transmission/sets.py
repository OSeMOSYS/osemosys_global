"""Function to create sets."""

from data import format_transmission_name

def create_tech_set_trn(df_lines):

    df = df_lines[['From', 'To']].copy()
    df = format_transmission_name(df)
    df = df.rename(columns = {'TECHNOLOGY' : 'VALUE'})    
    
    return df