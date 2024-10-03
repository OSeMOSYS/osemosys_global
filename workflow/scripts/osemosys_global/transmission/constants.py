"""Constants for the transmission module"""

"""Set iar and oar values for custom transmission entries. I.e. oar of 0.9
assumes 10% losses. Relevant for the user_defined_capacity function."""
DF_IAR_CUSTOM_VAL = 1
DF_OAR_CUSTOM_VAL = 0.9

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