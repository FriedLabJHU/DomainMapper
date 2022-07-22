import os
import sys
import shutil
from xml import dom
from setuptools import setup, find_packages

if os.path.exists('build') == True:
	print("build exists")
	shutil.rmtree('./build')

try:
    import requests
except ImportError:
    sys.stderr.write('requests is not installed, you can find it at: https://docs.python-requests.org/en/latest/ \n')
    sys.exit()

if [int(dgt) for dgt in requests.__version__.split('.')[:2]] < [2, 0]:
    sys.stderr.write('DomainMapper requires requests v2.0 or later, you can find it at: https://docs.python-requests.org/en/latest/ \n')
    sys.exit()

try:
    import scipy
except ImportError:
    sys.stderr.write('SciPy is not installed, you can find it at: https://scipy.org/ \n')
    sys.exit()

if [int(dgt) for dgt in scipy.__version__.split('.')[:2]] < [1, 6]:
    sys.stderr.write('DomainMapper requires scipy v1.6 or later, you can find it at: https://scipy.org/ \n')
    sys.exit()

try:
    import Bio
except ImportError:
    sys.stderr.write('BioPython is not installed, you can find it at: https://biopython.org/ \n')
    sys.exit()

if [int(dgt) for dgt in Bio.__version__.split('.')[:2]] < [1, 6]:
    sys.stderr.write('DomainMapper requires BioPython v1.6 or later, you can find it at: https://biopython.org/ \n')
    sys.exit()

try:
    import numpy
except ImportError:
    sys.stderr.write('numpy is not installed, you can find it at: https://numpy.org/ \n')
    sys.exit()

if [int(dgt) for dgt in numpy.__version__.split('.')[:2]] < [1, 17]:
    sys.stderr.write('DomainMapper requires numpy v1.17 or later, you can find it at: https://numpy.org/ \n')
    sys.exit()

if sys.version_info[:2] < (3, 4):
    print("DomainMapper requires Python 3.4 or later. Python {}.{} detected".format(*sys.version_info[:2]))
    print("Please upgrade your version of Python.")
    sys.exit(-1)

setup(
    name = "DomainMapper",
    version = "3.0.0",
    author = "Edgar Manriquez-Sandoval - Fried Lab - JHU",
    author_email = "emanriq1@jhu.edu",
    url="https://github.com/FriedLabJHU/DomainMapper",
    description = ("A parser for hmmscan full outputs built around ECOD domain definitions"),
    packages=["DomainMapper"],
    package_dir = {'':'src'},
    install_requires=[
        "requests",
        "bio",
        "numpy",
        "scipy"
    ],
    package_data={'': ['ecod.latest.domains.npy']},
    include_package_data=True,
    scripts=['src/DomainMapper/dommap'],
    python_requires=">=3.5"
)