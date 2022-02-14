"""Runs tests against long and short OSeMOSYS implementations and various solvers

Checks::

- OSeMOSYS Short formulation
- OSeMOSYS Long formulation
- OSeMOSYS Fast formulation

and solvers::

- GLPSOL
- CBC [awaiting release of otoole v0.7]

"""
import os
from pytest import fixture
import pytest
import pandas as pd
from typing import Tuple, Any

import tempfile

from subprocess import run


@fixture()
def get_folder():
    return os.path.dirname(os.path.abspath(__file__))


@fixture(
    scope="function",
    params=[("glpk", "long"), ("glpk", "short"), (("glpk", "fast"))],
    ids=["glpk-long", "glpk-short", "glpk-fast"],
)
def run_model(request, tmpdir, get_folder) -> Tuple[Any, str]:
    """This parameterised pytest fixture is used to run the different OSeMOSYS formulations

    The short, long and fast versions of the code are run for GLPK

    With the addition of a post-processing script for CBC, it will be possible to include
    CBC in the same tests.

    """

    if request.param[1] == "long":
        model_file = os.path.join(get_folder, "../src/osemosys.txt")
    elif request.param[1] == "short":
        model_file = os.path.join(get_folder, "../src/osemosys_short.txt")
    elif request.param[1] == "fast":
        model_file = os.path.join(get_folder, "../src/osemosys_fast.txt")

    results_folder = str(tmpdir.mkdir("results"))

    os.chdir(os.path.join(results_folder, ".."))

    data_file = os.path.join(get_folder, "utopia.txt")

    if request.param[0] == "glpk":

        arguments = ["glpsol", "-m", model_file, "-d", data_file]
        output = run(arguments, capture_output=True, text=True)
        return output, os.path.join(results_folder)

    elif request.param[0] == "cbc":

        lp_file = tempfile.NamedTemporaryFile(suffix=".lp").name

        arguments = [
            "glpsol",
            "-m",
            model_file,
            "-d",
            data_file,
            "--wlp",
            lp_file,
            "--check",
        ]
        run(arguments)

        cbc_results_file = os.path.join(results_folder, "results.txt")

        arguments = ["cbc", lp_file, "-sec", "15", "solve", "-solu", cbc_results_file]
        output = run(arguments)
        # Add post-processing for CBC results here
        return output, cbc_results_file


class TestOsemosysOutputs:
    def test_mathprog_run_normal(self, run_model):
        output, results_folder = run_model
        try:
            assert "OPTIMAL LP SOLUTION FOUND" in output.stdout
        except AssertionError as ex:
            raise AssertionError(str(ex), output)
        assert "obj =   2.944686269e+04" in output.stdout

    def test_results_exist(self, run_model):
        """OSeMOSYS short and long versions should produce 29 CSV files

        These are placed into the folder specified by the ``ResultPath``
        parameter
        """
        _, results_folder = run_model
        print(results_folder)
        assert os.path.exists(str(results_folder))
        print(os.listdir(str(results_folder)))
        assert len(os.listdir(str(results_folder))) == 30

    def test_results_read_accumulated_new_capacity(self, run_model):
        """
        """
        _, results_folder = run_model
        result_file = os.path.join(str(results_folder), "AccumulatedNewCapacity.csv")
        actual = pd.read_csv(result_file)
        actual = actual[actual.YEAR == 2010].reset_index(drop=True)

        expected = pd.DataFrame(
            columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"],
            data=[
                ["UTOPIA", "E01", 2010, 2.279801],
                ["UTOPIA", "E31", 2010, 0.110000],
                ["UTOPIA", "IMPDSL1", 2010, 77.597496],
                ["UTOPIA", "IMPHCO1", 2010, 191.565506],
                ["UTOPIA", "RHE", 2010, 46.867723],
                ["UTOPIA", "RHO", 2010, 46.135248],
                ["UTOPIA", "RL1", 2010, 18.901890],
                ["UTOPIA", "SRE", 2010, 0.100000],
                ["UTOPIA", "TXD", 2010, 11.690000],
                ["UTOPIA", "RIV", 2010, 5.587785],
            ],
        )

        pd.testing.assert_frame_equal(actual, expected)

    def test_results_read_new_capacity(self, run_model):
        """
        """
        _, results_folder = run_model

        result_file = os.path.join(str(results_folder), "NewCapacity.csv")
        actual = pd.read_csv(result_file)
        actual = (
            actual.groupby(by=["REGION", "TECHNOLOGY"], as_index=False)
            .sum()
            .drop(columns="YEAR")
        )

        expected = pd.DataFrame(
            columns=["REGION", "TECHNOLOGY", "VALUE"],
            data=[
                ["UTOPIA", "E01", 2.279801],
                ["UTOPIA", "E31", 0.110000],
                ["UTOPIA", "IMPDSL1", 1717.805326],
                ["UTOPIA", "IMPHCO1", 1451.847478],
                ["UTOPIA", "RHE", 46.867723],
                ["UTOPIA", "RHO", 46.135248],
                ["UTOPIA", "RIV", 97.919280],
                ["UTOPIA", "RL1", 34.303990],
                ["UTOPIA", "SRE", 0.100000],
                ["UTOPIA", "TXD", 17.790000],
            ],
        )

        pd.testing.assert_frame_equal(actual, expected)

    def test_results_read_investment_cost(self, run_model):
        """
        """
        _, results_folder = run_model

        result_file = os.path.join(str(results_folder), "CapitalInvestment.csv")
        df = pd.read_csv(result_file)

        actual = (
            df.groupby(by=["REGION", "TECHNOLOGY"], as_index=False)
            .sum()
            .drop(columns="YEAR")
        )

        expected = pd.DataFrame(
            columns=["REGION", "TECHNOLOGY", "VALUE"],
            data=[
                ["UTOPIA", "E01", 2845.536735],
                ["UTOPIA", "E31", 330.000000],
                ["UTOPIA", "RHE", 4218.095050],
                ["UTOPIA", "RHO", 4613.524752],
                ["UTOPIA", "SRE", 10.000000],
                ["UTOPIA", "TXD", 18572.760000],
            ],
        )

        pd.testing.assert_frame_equal(actual, expected)

    def test_results_read_varom(self, run_model):
        """
        """
        _, results_folder = run_model

        result_file = os.path.join(
            str(results_folder), "AnnualVariableOperatingCost.csv"
        )
        df = pd.read_csv(result_file)

        actual = (
            df.groupby(by=["REGION", "TECHNOLOGY"], as_index=False)
            .sum()
            .drop(columns="YEAR")
        )

        expected = pd.DataFrame(
            columns=["REGION", "TECHNOLOGY", "VALUE"],
            data=[
                ["UTOPIA", "E01", 90.471708],
                ["UTOPIA", "E31", 0.000310],
                ["UTOPIA", "E51", 0.000563],
                ["UTOPIA", "IMPDSL1", 10618.931276],
                ["UTOPIA", "IMPHCO1", 1884.827246],
                ["UTOPIA", "RHE", 0.002048],
                ["UTOPIA", "RHO", 0.006237],
                ["UTOPIA", "RIV", 0.000968],
                ["UTOPIA", "RL1", 0.001841],
                ["UTOPIA", "TXD", 0.001709],
            ],
        )

        pd.testing.assert_frame_equal(
            actual, expected, check_exact=False, rtol=1e-3
        )

    def test_results_emissions(self, run_model):
        """
        """
        _, results_folder = run_model

        result_file = os.path.join(str(results_folder), "AnnualEmissions.csv")
        df = pd.read_csv(result_file)

        actual = (
            df.groupby(by=["REGION", "EMISSION"], as_index=False)
            .sum()
            .drop(columns="YEAR")
        )

        expected = pd.DataFrame(
            columns=["REGION", "EMISSION", "VALUE"],
            data=[["UTOPIA", "CO2", 163.516797], ["UTOPIA", "NOX", 170.895000]],
        )

        pd.testing.assert_frame_equal(actual, expected)

    def test_total_tech_model_period_activity(self, run_model):
        """
        """

        _, results_folder = run_model

        result_file = os.path.join(str(results_folder), "TotalTechnologyModelPeriodActivity.csv")
        actual = pd.read_csv(result_file)

        columns = ["REGION", "TECHNOLOGY", "VALUE"]
        data = [
            ["UTOPIA", "E01", 301.572359341649],
            ["UTOPIA", "E31", 30.98719154944],
            ["UTOPIA", "E51", 56.29176],
            ["UTOPIA", "IMPDSL1", 1061.89312756574],
            ["UTOPIA", "IMPHCO1", 942.413622942653],
            ["UTOPIA", "RHE", 204.751310891089],
            ["UTOPIA", "RHO", 623.698689108911],
            ["UTOPIA", "RL1", 184.1],
            ["UTOPIA", "TXD", 170.895],
            ["UTOPIA", "RIV", 96.8349735920001],
            ]

        expected = pd.DataFrame(
            columns=columns,
            data=data,
        )

        pd.testing.assert_frame_equal(actual, expected)

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_results_rate_of_production(self, run_model):
        """REGION,TIMESLICE,TECHNOLOGY,FUEL,YEAR,VALUE
        """
        _, results_folder = run_model

        result_file = os.path.join(str(results_folder), "RateOfProductionByTechnology.csv")
        df = pd.read_csv(result_file)

        actual = (
            df.groupby(by=["REGION", "TIMESLICE", "TECHNOLOGY", "FUEL"], as_index=False)
            .sum()
            .drop(columns="YEAR")
        )

        expected = pd.DataFrame(
            columns=["REGION", "TIMESLICE", "TECHNOLOGY", "FUEL", "VALUE"],
            data=[["UTOPIA", "WS", "TECH", "FUEL", 163.516797]],
        )

        pd.testing.assert_frame_equal(actual, expected)
