"""
Test script for pywopwop - https://github.com/fchirono/pywopwop

    This script creates a toy dataset manually, containing a

    rites this data to a set
    of geometry and loading PSU-WOPWOP v1.0 files.


Author:
    Fabio Casagrande Hirono
    Oct 2021
"""

import numpy as np
import pywopwop as PWW

# %%
# create toy dataset - random numbers for mesh coordinates, normal vectors and
# surface pressures

name = 'Zone 1'

iMax = 20
jMax = 61


rng = np.random.default_rng()

XYZ_coord = rng.normal(size=(3, iMax, jMax))
normal_coord = rng.normal(size=(3, iMax, jMax))
loading_data = rng.normal(size=(iMax, jMax))


# %%
# initialize new instance
myWopwopData = PWW.PWWPatch()

myWopwopData.geom_type = 'geometry'
myWopwopData.is_structured = True
myWopwopData.geometry_time_type = 'constant'
myWopwopData.centered_type = 'node'
myWopwopData.float_type = 'single'
myWopwopData.iblank_included = 'n'
myWopwopData.units_string = 'Pa'

# 1024-byte comment string for geometry file
myWopwopData.geometry_comment = 'Test file for pywopwop - geometry'


myWopwopData.loading_time_type = 'constant'
myWopwopData.loading_data_type = 'surf_pressure'
myWopwopData.loading_ref_frame = 'blade_fixed'
# 1024-byte comment string for loading file
myWopwopData.loading_comment = 'Test file for pywopwop - loading'


# adds zone 1 info - mesh coords, normal coords, loading data
myWopwopData.add_StructuredZone(name, XYZ_coord, normal_coord,
                                calc_thickness_noise=True,
                                loading_data=loading_data)

# add zone 2 - no loading data
myWopwopData.add_StructuredZone(name[:-1] + '2', XYZ_coord+1, normal_coord+1,
                                calc_thickness_noise=True)

# add zone 3 - loading data, no thickness noise
myWopwopData.add_StructuredZone(name[:-1] + '3', XYZ_coord+2, normal_coord+2,
                                calc_thickness_noise=False,
                                loading_data=loading_data+2)

myWopwopData.print_info()

# %% write files

myWopwopData.write_geometry_file('BladeGeometry.dat')
myWopwopData.write_loading_file('BladeLoading.dat')
