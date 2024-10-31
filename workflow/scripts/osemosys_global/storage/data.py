"""Functions to extract and format relevent data for storage."""

def get_years(start: int, end: int) -> range:
    return range(start, end + 1)

def get_years_build_rates(row):
    row['YEAR'] = get_years(row['START_YEAR'], row['END_YEAR'])
    
    return row

