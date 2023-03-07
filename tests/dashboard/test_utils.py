"""Module for testing dashboard utils"""

import osemosys_global.dashboard.utils as utils
import pandas as pd
from pytest import mark, fixture
from pandas.testing import assert_frame_equal


class TestAddDefaultTemporalValue:
    
    def years():
        return [2020, 2021, 2022, 2023, 2024, 2025]
    
    def timeslices():
        return ["S1D1", "S1D2", "S2D1", "S2D2", "S3D1", "S3D2"]
    
    def default_value():
        return 0
    
    df_year_in = pd.DataFrame([
        ["TECH", 2020, 100],
        ["TECH", 2021, 200],
        ["TECH", 2023, 300],
        ["TECH", 2025, 400],
    ], columns=["TECHNOLOGY", "YEAR", "VALUE"])
    
    df_ts_in = pd.DataFrame([
        ["TECH", "S1D1", 2020, 100],
        ["TECH", "S2D1", 2020, 200],
        ["TECH", "S2D2", 2020, 300],
        ["TECH", "S3D2", 2020, 400],
    ], columns=["TECHNOLOGY", "TIMESLICE", "YEAR", "VALUE"])
    
    df_year_out = pd.DataFrame([
        ["TECH", 2020, 100],
        ["TECH", 2021, 200],
        ["TECH", 2022, 0],
        ["TECH", 2023, 300],
        ["TECH", 2024, 0],
        ["TECH", 2025, 400],
    ], columns=["TECHNOLOGY", "YEAR", "VALUE"])

    df_ts_out = pd.DataFrame([
        ["TECH", "S1D1", 2020, 100],
        ["TECH", "S1D2", 2020, 0],
        ["TECH", "S2D1", 2020, 200],
        ["TECH", "S2D2", 2020, 300],
        ["TECH", "S3D1", 2020, 0],
        ["TECH", "S3D2", 2020, 400],
    ], columns=["TECHNOLOGY", "TIMESLICE", "YEAR", "VALUE"])
    
    test_data = [
        (df_year_in, df_year_out, "YEAR", years(), default_value()),
        (df_ts_in, df_ts_out, "TIMESLICE", timeslices(), default_value())
    ]
    
    @mark.parametrize("df, expected, column, indices, default_value", test_data, ids=["year", "timeslice"])
    def test_add_default_temporal_value(self, df, expected, column, indices, default_value):
        actual = utils.add_default_values(
            df=df,
            column=column,
            default_indices=indices,
            default_value=default_value
        )
        assert_frame_equal(actual, expected)
        
        