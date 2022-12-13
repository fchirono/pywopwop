"""
pywopwop - https://github.com/fchirono/pywopwop

    This script reads a pair of geometry and loading PSU-WOPWOP files, and
    rewrites the same data using different file names.

    The files can then be compared for bitwise equality in bash using:

        cmp /path/to/filename1 /path/to/filename2


Author:
    Fabio Casagrande Hirono
    Nov 2022
"""

import numpy as np
rng = np.random.default_rng()

import pywopwop as PWW

write_files = True

# %% Define directory containing geometry and loading files, plus filenames

# *****************************************************************************
# # CASE 1: constant geometry blade, constant compact loading vectors [OK]
# wopwop_dir = '../../OneDrive - University of Salford/PSU-WOPWOP_v3.4.4/case1/'
# geometry_filename = 'gyrodyne.dat'
# loading_filename = 'gyrodyneLoading.dat'


# *****************************************************************************
# # CASE 2: permeable surface with constant geometry, periodic flow data [OK]
# #   --> periodic loading does NOT have extra time step at the end [?]
# wopwop_dir = '../../OneDrive - University of Salford/PSU-WOPWOP_v3.4.4/case2/'
# geometry_filename = 'reduced_geometry.dat'
# loading_filename = 'reduced_flow_data.dat'


# *****************************************************************************
# CASE 4: constant surface geometry, periodic loading vectors [OK]
#   --> periodic loading has an extra time step at the end [?]
wopwop_dir = '../../OneDrive - University of Salford/PSU-WOPWOP_v3.4.4/case4/'
geometry_filename = 'HART_Blade.dat'
loading_filename = 'HART_Bladeloading.dat'


# *****************************************************************************
# # CASE 5: constant geometry, aperiodic loading [OK]
# wopwop_dir = '../../OneDrive - University of Salford/PSU-WOPWOP_v3.4.4/case5/'
# geometry_filename = 'constGeo_short.dat'
# loading_filename = 'AperLoadingShort.dat'


# %% initialize new instance, read files
myWopwopData = PWW.PWWPatch()
myWopwopData.read_geometry_file(wopwop_dir + geometry_filename)
myWopwopData.read_loading_file(wopwop_dir + loading_filename)

myWopwopData.print_info(zones_info=True)

# %% rewrite same data using different filenames
if write_files:
    myWopwopData.write_geometry_file(wopwop_dir + 'geometry_file.dat')
    myWopwopData.write_loading_file(wopwop_dir + 'loading_file.dat')

# re-read files
myWopwopData2 = PWW.PWWPatch()
myWopwopData2.read_geometry_file(wopwop_dir + 'geometry_file.dat')
myWopwopData2.read_loading_file(wopwop_dir + 'loading_file.dat')
myWopwopData2.print_info(zones_info=True)

# compare patches' contents
PWW.compare_pwwpatches(myWopwopData, myWopwopData2)
