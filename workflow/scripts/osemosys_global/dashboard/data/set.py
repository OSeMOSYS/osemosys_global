from __future__ import annotations
from dataclasses import dataclass

import pandas as pd

@dataclass
class SetData:
    """Holds set data"""
    
    _data = pd.DataFrame
    
    @property
    def years(self) -> list[int]:
        return self._data["VALUE"].to_list()
    
    @property
    def country_codes(self) -> list[str]:
        pass