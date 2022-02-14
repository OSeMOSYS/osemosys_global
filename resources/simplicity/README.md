# OSeMOSYS Example Model: Simplicity

[![DOI](https://zenodo.org/badge/214192147.svg)](https://zenodo.org/badge/latestdoi/214192147)
[![goodtables.io](https://goodtables.io/badge/github/OSeMOSYS/simplicity.svg)](https://goodtables.io/github/OSeMOSYS/simplicity)

This is an example OSeMOSYS model published as a Frictionless Data [Tabular Data Package](https://frictionlessdata.io/specs/tabular-data-package/).

The Data Package was constructed with the help of the Python **otoole** package available on [PyPI](https://pypi.org/project/otoole/) and [Github](https://github.com/OSeMOSYS/otoole).

Simplicity `v0.2` requires `v0.5.4` or later of **otoole**.
You can use **otoole** to generate a GNU MathProg data file from the dataset with the following commands and then run OSeMOSYS.

```bash
# Install the OSeMOSYS toolkit
pip install otoole>=0.5.4
# Download the dataset and build a GNU MathProg datafile 
#   On Linux or OSX use wget
wget https://zenodo.org/record/3707794/files/OSeMOSYS/simplicity-v0.2.1.zip
unzip simplicity-v0.2.1.zip -d simplicity
#   Otherwise download manually and unzip, and then...
otoole convert datapackage datafile simplicity/OSeMOSYS-simplicity-11a3a26/datapackage.json ./simplicity.txt
# Solve the model
glpsol -m OSeMOSYS.txt -d simplicity.txt
```
