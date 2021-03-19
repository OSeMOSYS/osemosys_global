
# -*- coding: utf-8 -*-
from pkg_resources import DistributionNotFound, get_distribution
from OPG_geographic_filter import extract_country

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

__all__ = ["extract_country"]

