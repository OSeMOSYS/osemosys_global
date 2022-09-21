## config.yaml

This configuration file contains all parameters that the user can edit. The table below describes the available parameters and their limits. 

| Parameter  | Description | Limits |
|------------|-------------|--------|
| `scenario` | Scenario name |     |
| `startYear` | Start year of model | 2015 |
| `endYear`| End year of model | 2100 |
| `dayType`| To be implemented  |   |
| `dayParts` | Hours included in each day part | Include all values between 0-24 |
| `timeshift` | Shifts location of modelled time zero relative to UTC | Between -11 and 12 |
| `seasons` | Months included in each season | Include all values between 1-12 |
| `geographic_scope` | Countries to include in scenario | Three letter codes using [this source](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) |
| `crossborderTrade` | Enable trading of resources between all countries in scenario | True or False |
| `emission_penalty`   | Global carbon tax in Million Dollars per MegaTonne | Number greater than or equal to zero |
| `no_invest_technologies` | Technologies that can not be invested in | Technology codes from list below |
| `results_by_country` | Plot results by country in addition to system level results | True or False |
| `solver` | LP Solver to use | `cbc`, `gurobi`, `cplex` |
| `user_defined_capacity` | Modelled capacity additions |  |

