"""Constants for the demand module"""

# Sets the spatial resolution for the regession, right now can only be applied to continental level by setting 'child_object'.
# Note that not all countries have historical data available so country-level regression can't be applied for all.
SPATIAL_RESOLUTION = "child_object"

# Include urbanization as part of the regression? 'Yes' or 'No'
URBANIZATION = "Yes"

# Set which SSP data and sources are to be used for the country-level GDP|PPP and Population projects.
PATHWAY = "SSP2"
POP_COUNTRIES_SOURCE = "IIASA-WiC POP"  # Options are 'IIASA-WiC POP' and 'NCAR'
GDP_PPP_COUNTRIES_SOURCE = "OECD Env-Growth"  # Options are 'IIASA GDP' and 'OECD Env-Growth'. 'OECD Env-Growth' has more entries.
URB_COUNTRIES_SOURCE = "NCAR"  #'NCAR' is the only option.

# Projection range
START_YEAR = 2010
END_YEAR = 2100
