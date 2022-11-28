"""
Test script for pywopwop - https://github.com/fchirono/pywopwop

    This script creates a toy dataset manually and shows how to introduce
    data into a PWWPatch instance.

    The example dataset is a collection of structured meshes, with constant
    geometry and constant loading. It contains random values for mesh
    coordinates, normal vectors and surface pressures.


Author:
    Fabio Casagrande Hirono
    Nov 2022
"""

import numpy as np
import pywopwop as PWW

write_files = False

# %% create toy dataset with structured geometry

name = 'Zone 1'
iMax = 20
jMax = 61

# uses random numbers for mesh coordinates, normal vectors and surface pressures
rng = np.random.default_rng()

XYZ_coord = rng.normal(size=(3, iMax, jMax))
normal_coord = rng.normal(size=(3, iMax, jMax))
loading_data = rng.normal(size=(iMax, jMax))


# %% initialize PWWPatch instance
myWopwopData = PWW.PWWPatch()

# add general header info
myWopwopData.is_structured = True
myWopwopData.centered_type = 'node'
myWopwopData.float_type = 'single'
myWopwopData.has_iblank = False
myWopwopData.set_units_string('Pa')

# add geometry info
myWopwopData.geometry_type = 'geometry'
myWopwopData.geometry_time_type = 'constant'
myWopwopData.set_geometry_comment('Test file for pywopwop - geometry')

# add loading info
myWopwopData.loading_time_type = 'constant'
myWopwopData.loading_data_type = 'surf_pressure'
myWopwopData.loading_ref_frame = 'blade_fixed'
myWopwopData.set_loading_comment('Test file for pywopwop - loading')


# %% Adds structured zones, one at a time

# *****************************************************************************
# adds zone 1 info - name: 'Zone 1'
#   --> adds random mesh coords, normal coords, loading data
#   --> has thickness noise
#   --> has loading data

myWopwopData.add_StructuredZone(name,
                                XYZ_coord, normal_coord,
                                calc_thickness_noise=True,
                                loading_data=loading_data)

# *****************************************************************************
# adds zone 2 info - name: 'Zone 2'
#   --> XYZ coords and normal coords are identical to zone 1, plus 1 unit
#   --> has thickness noise
#   --> doesn't have loading data

myWopwopData.add_StructuredZone(name[:-1] + '2',
                                XYZ_coord+1, normal_coord+1,
                                calc_thickness_noise=True)

# *****************************************************************************
# adds zone 1 info - name: 'Zone 3'
#   --> XYZ coords, normal coords, and loading data are identical to zone 1,
#       plus 2 units
#   --> doesn't have thickness noise
#   --> has loading data

myWopwopData.add_StructuredZone(name[:-1] + '3',
                                XYZ_coord+2, normal_coord+2,
                                calc_thickness_noise=False,
                                loading_data=loading_data+2)

# %% # print info, write files

# print info about PWWPatch, including zones' info
myWopwopData.print_info(zones_info=True)

if write_files:
    myWopwopData.write_geometry_file('BladeGeometry.dat')
    myWopwopData.write_loading_file('BladeLoading.dat')
