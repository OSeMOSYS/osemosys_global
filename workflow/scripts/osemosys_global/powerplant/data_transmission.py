"""Functions to extract and format relevent data for tranmission."""

def format_transmission_name(df):
    '''Formats PLEXOS transmission names into OSeMOSYS Global names.

    Args:
        :param df: Pandas DataFrame with columns 'From' and 'To' describing the 
               transmission from and to contries. ie. 
    
    Returns: 
        :param df: Same as df_in, except the 'From' and 'To' columns are replaced 
            with a single 'TECHNOLOGY' column holding OSeMOSYS Global 
            naming conventions 

    Example:
        df = pd.DataFrame(
            [[AF-COD, AF-COG, 0.001],
            [EU-AUT, EU-SVK, 0.004],
            [AS-LBN, AS-SYR, 0.006]], 
            columns = ['From', 'To', 'Losses']
        )
        pd.DataFrame(
            [[0.001,TRNCODXXCOGXX],
            [0.004,TRNAUTXXSVKXX],
            [0.006,TRNLBNXXSYRXX]] 
            columns = ['Losses','TECHNOLOGY'])'''

    # If from column has length 6 then it's the last three chars plus XX
    df.loc[df["From"].str.len() == 6, "From"] = (df["From"].str[3:6] + "XX")

    # If from column has length 9 then it's the 3:6 and 7:9 three chars plus XX
    df.loc[df["From"].str.len() == 9, "From"] = (
        df["From"].str[3:6] + df["From"].str[7:9])

    # If to column has length 6 then it's the last three chars plus XX
    df.loc[df["To"].str.len() == 6, "To"] = (df["To"].str[3:6] + "XX")

    # If to column has length 9 then it's the 3:6 and 7:9 three chars plus XX
    df.loc[df["To"].str.len() == 9, "To"] = (
        df["To"].str[3:6] + df["To"].str[7:9])

    # Combine From and To columns.
    df["TECHNOLOGY"] = ("TRN" + df["From"] + df["To"])

    # Drop to and from columns
    df = df.drop(["From", "To"], axis=1)

    return df