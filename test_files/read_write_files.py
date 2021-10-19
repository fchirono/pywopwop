"""
Test script for pywopwop - https://github.com/fchirono/pywopwop

    This script reads a pair of geometry and loading PSU-WOPWOP files, and
    rewrites the same data using different file names.

    The files can then be compared for bitwise equality in bash using:

        cmp /path/to/filename1 /path/to/filename2


Author:
    Fabio Casagrande Hirono
    Oct 2021
"""

import numpy as np
rng = np.random.default_rng()

import pywopwop as PWW


# %% Define directory containing geometry and loading files, plus filenames

# wopwop_dir = '../../../PSU-WOPWOP_v3.4.4/case1/'
# geometry_filename = 'gyrodyne.dat'
# loading_filename = 'gyrodyneLoading.dat'

# aperiodic loading not supported yet - loading file will not be identical
wopwop_dir = '../../../PSU-WOPWOP_v3.4.4/case5/'
geometry_filename = 'constGeo_short.dat'
loading_filename = 'AperLoadingShort.dat'

# %%
# initialize new instance, read files
myWopwopData = PWW.PWWPatch()
myWopwopData.read_geometry_file(wopwop_dir + geometry_filename)
myWopwopData.read_loading_file(wopwop_dir + loading_filename)

myWopwopData.print_info()

# rewrite same data in different filenames
myWopwopData.write_geometry_file('geometry_file.dat')
myWopwopData.write_loading_file('loading_file.dat')