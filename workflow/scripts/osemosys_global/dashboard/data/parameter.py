from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

import pandas as pd

@dataclass
class SourceData:
    """Holds parameter data"""
    
    _data = pd.DataFrame
    
    def filter_for_bar():
        pass
    
    def filer_for_pie():
        pass
    
    def filter_for_line():
        pass
    
    @property
    def row_count(self) -> int:
        return self._data.shape[0]

