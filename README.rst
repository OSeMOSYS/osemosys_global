===============
osemosys_global
===============

A global power system model generator for OSeMOSYS

Installation
------------

Due to the data distributed as part of the package, an installation through PyPi is not
available. Instead, please install using the following instructions.



python setup.py develop


Usage
-----


Usage::

    import osemosys_global as og
    og.get_data()
    og.extract_country()