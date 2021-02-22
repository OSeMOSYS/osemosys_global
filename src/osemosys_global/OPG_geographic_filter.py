#!/usr/bin/env python
# coding: utf-8

"""Filter osemosys_global datapackaged based on user-defined geographic scope

"""
import pandas as pd
import os
import sys
from typing import List


def extract_country(iso_country_code: List[str], model_name: str):
    """Filter one or more countries from the global data set

    Arguments
    ---------
    iso_country_code: List[str]
        List of three-letter ISO country codes
    model_name: str
        Folder name into which to insert data
    """
    if not os.path.exists(os.path.join(os.getcwd(),
                            model_name,
                            'data')):
        os.makedirs(os.path.join(os.getcwd(),
                            model_name,
                            'data'))

    for each_csv in (os.listdir(os.path.join(os.getcwd(),
                                             'osemosys_global_model', 'data'))):

        df = pd.read_csv(os.path.join(os.getcwd(),
                                      'osemosys_global_model', 'data',
                                      each_csv)
                        )
        if not df.empty:
            if 'TECHNOLOGY' in df.columns:
                df = df.loc[df['TECHNOLOGY'].str[3:6].isin(iso_country_code) |
                            df['TECHNOLOGY'].str[6:9].isin(iso_country_code) |
                            df['TECHNOLOGY'].str[8:11].isin(iso_country_code)]

            if 'FUEL' in df.columns:
                df = df.loc[df['FUEL'].str[3:6].isin(iso_country_code) |
                            df['FUEL'].str[6:9].isin(iso_country_code)]

            if each_csv == 'FUEL.csv':
                df = df.loc[df['VALUE'].str[3:6].isin(iso_country_code) |
                            df['VALUE'].str[6:9].isin(iso_country_code)]

            if each_csv == 'TECHNOLOGY.csv':
                df = df.loc[df['VALUE'].str[3:6].isin(iso_country_code) |
                            df['VALUE'].str[6:9].isin(iso_country_code) |
                            df['VALUE'].str[8:11].isin(iso_country_code)]


        df.to_csv(os.path.join(os.getcwd(),
                               model_name,
                               'data',
                               each_csv),
                  index=None)


if __name__ == '__main__':

    args = sys.argv[1:]
    if len(args) != 2:
        print("Usage: python OPG_geographic_filter <country_code> <output_name>")
        exit(1)

    iso_country_code = args[0]
    model_name = args[1]

    extract_country([iso_country_code], model_name)