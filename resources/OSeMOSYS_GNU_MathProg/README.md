# OSeMOSYS GNU MathProg

Thanks for using OSeMOSYS and welcome to the OSeMOSYS community.

To run OSeMOSYS, enter the following line into your command prompt and
data file name:

    glpsol -m osemosys.txt -d  ../Training_Case_Studies/utopia.txt -o results.csv

Alternatively, install GUSEK (http://gusek.sourceforge.net/gusek.html)
and run the model within this integrated development environment (IDE).
To do so, open the datafile (e.g. `utopia.txt`) and
select "Use External .dat file" from the Options menu.
Then change to the model file and select the "Go" icon or press F5.

## Developers - Testing

This repository uses Travis CI to run regression tests and
harmonisation tests across each of the OSeMOSYS GNU MathProg normal and short
implementations.

Each push to a branch on the repository, or submission of a pull
request triggers a build on Travis CI, with the corresponding status reported
back in the pull request comments.

The tests must pass before a pull request may be merged into the main
repository.

Tests are defined using the Python package ``pytest`` and the runs are
configured within the Travis CI configuration file ``.travis.yml``.

The tests are stored in the ``tests`` folder.

### Running the tests

To run the tests on your local computer, you need a Python 3.7 installation.
The easiest way to install this is using
[miniconda](https://docs.conda.io/en/latest/miniconda.html).

Then you need to install pytest `conda install pytest pandas` and can then run the tests
using the command `pytest`.

Each of the tests in the `tests` folder runs an OSeMOSYS model file and checks that
the output matches a given value.

## Developers - Creating a new release and deploying to Github Releases

Creating a new release for OSeMOSYS GNU MathProg is as simple as creating a new semver
compliant tag and pushing the tag to a branch.  Travis CI will then run the tests, and 
if they pass, create the package using the `makefile` found in the root of the repository.
The makefile includes the contents of the `src` and `scripts` folders, an html render of 
`src/README.md` and then zips them up deploying them to Github Releases and providing 
the contents of `docs/changelog.md` as a release description.

### 1. Update the changelog

Add a new heading with the version number you will use and include a description 
of the changes since the last release.

It can be useful to view the log of git commits since the previous release using the 
following command:

    git log <yourlasttag>..HEAD

### 2. Run the tests locally

Follow the instructions provided to run the tests.

### 3. Create a new tag and push to Github

Create a new annotated tag:

    git tag -a v1.0.0 -m "A descriptive message for the release"

Please follow the Semantic Versioning [guidelines](https://semver.org/).
To create an alpha or beta release (pre-release) you would do the following:

    git tag -a v1.0.0-alpha.1 -m "An alpha release"
    git tag -a v1.0.0-beta.1 -m "An beta release"
    git tag -a v1.0.0 -m "First stable official release!"

### 4. Check that the package is deployed successfully

You can follow the release process at the 
[Travis CI service](https://travis-ci.com/github/OSeMOSYS/OSeMOSYS_GNU_MathProg/branches).

Finally, check that the release appears on the [Github Releases page](https://github.com/OSeMOSYS/OSeMOSYS_GNU_MathProg/releases).
