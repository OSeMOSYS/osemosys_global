"""Constants for plotting"""

MONTH_NAMES = {
    1: 'Jan',
    2: 'Feb',
    3: 'Mar',
    4: 'Apr',
    5: 'May',
    6: 'Jun',
    7: 'Jul',
    8: 'Aug',
    9: 'Sep',
    10: 'Oct',
    11: 'Nov',
    12: 'Dec',
}

DAYS_PER_MONTH = {
    'Jan': 31,
    'Feb': 28,
    'Mar': 31,
    'Apr': 30,
    'May': 31,
    'Jun': 30,
    'Jul': 31,
    'Aug': 31,
    'Sep': 30,
    'Oct': 31,
    'Nov': 30,
    'Dec': 31,
}

COLORS = {
    "BIO":"darkgreen",
    "CCG":"lightcoral",
    "COA":"black",
    "COG":"peru",
    "CSP":"wheat",
    "ELC":"gold",
    "GAS":"orange",
    "GEO":"darkseagreen",
    "HYD":"dodgerblue",
    "OCG":"firebrick",
    "OIL":"lightgrey",
    "OTH":"teal",
    "PET":"grey",
    "SOL":"gold",
    "SPV":"gold",
    "URN":"mediumseagreen",
    "WAS":"darkkhaki",
    "WAV":"navy",
    "WOF":"violet",
    "WON":"blueviolet",
    "INT":"darkgreen",
}

SET_DTYPES = {
    "DAILYTIMEBRACKET": int,
    "EMISSION":str,
    "FUEL":str,
    "MODE_OF_OPERATION":int,
    "REGION":str,
    "SEASON":str,
    "STORAGE":str,
    "TECHNOLOGY":str,
    "TIMESLICE":str,
    "YEAR":int,
}