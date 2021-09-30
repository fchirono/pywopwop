# -*- coding: utf-8 -*-
"""
Sample code for generating an airfoil surface, extracting the normal vectors
and visualising it using PyVista.

Based on the following issue: https://github.com/pyvista/pyvista-support/issues/506

Author:
    Fabio Casagrande Hirono
    Sep 2021
"""


import numpy as np
import pyvista as pv

import pywopwop as PWW

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# build a generic airfoil using 'airfoils' package (available on PyPI)
# ---> https://airfoils.readthedocs.io/en/latest/index.html

# from airfoils import Airfoil
# foil = Airfoil.NACA4('4812')
# # foil = Airfoil.NACA4('0012')

# Nx = 102
# x_foil = np.linspace(0, 1, Nx//2)
# y_foil_up = foil.y_upper(x_foil)
# y_foil_lo = foil.y_lower(x_foil)


# # concatenate top and bottom
# x = np.concatenate((x_foil, x_foil[::-1]))
# y = np.concatenate((y_foil_up, y_foil_lo[::-1]))

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Alternatively, define data points directly:

x = np.array([ 0.  , 0.02, 0.04, 0.06, 0.08, 0.1 , 0.12, 0.14, 0.16, 0.18, 0.2 ,
                0.22, 0.24, 0.26, 0.28, 0.3 , 0.32, 0.34, 0.36, 0.38, 0.4 , 0.42,
                0.44, 0.46, 0.48, 0.5 , 0.52, 0.54, 0.56, 0.58, 0.6 , 0.62, 0.64,
                0.66, 0.68, 0.7 , 0.72, 0.74, 0.76, 0.78, 0.8 , 0.82, 0.84, 0.86,
                0.88, 0.9 , 0.92, 0.94, 0.96, 0.98, 1.  , 1.  , 0.98, 0.96, 0.94,
                0.92, 0.9 , 0.88, 0.86, 0.84, 0.82, 0.8 , 0.78, 0.76, 0.74, 0.72,
                0.7 , 0.68, 0.66, 0.64, 0.62, 0.6 , 0.58, 0.56, 0.54, 0.52, 0.5 ,
                0.48, 0.46, 0.44, 0.42, 0.4 , 0.38, 0.36, 0.34, 0.32, 0.3 , 0.28,
                0.26, 0.24, 0.22, 0.2 , 0.18, 0.16, 0.14, 0.12, 0.1 , 0.08, 0.06,
                0.04, 0.02, 0.  ])

y = np.array([  0.        ,  0.03805971,  0.05349321,  0.06597081,  0.07664938,
                0.08600053,  0.09427433,  0.10162394,  0.10815206,  0.11393193,
                0.11901819,  0.12345299,  0.1272697 ,  0.13049523,  0.13315164,
                0.13525721,  0.13682722,  0.13787454,  0.13841001,  0.13844282,
                0.13797963,  0.13710874,  0.13595412,  0.13452314,  0.13282156,
                0.13085448,  0.12862644,  0.12614146,  0.12340309,  0.12041439,
                0.11717802,  0.11369624,  0.10997091,  0.10600351,  0.10179516,
                0.09734664,  0.09265833,  0.0877303 ,  0.08256226,  0.07715356,
                0.07150319,  0.06560979,  0.05947162,  0.05308659,  0.04645218,
                0.03956552,  0.03242328,  0.02502174,  0.0173567 ,  0.00942353,
                0.00121706, -0.0012992 ,  0.00117944,  0.00354851,  0.00580537,
                0.00794761,  0.00997303,  0.01187967,  0.01366579,  0.01532991,
                0.01687078,  0.01828745,  0.01957923,  0.0207457 ,  0.02178678,
                0.02270265,  0.02349386,  0.02416127,  0.02470611,  0.02512998,
                0.02543486,  0.02562315,  0.02569767,  0.02566172,  0.02551906,
                0.02527399,  0.02493133,  0.02449652,  0.02397561,  0.02337537,
                0.0227033 ,  0.0219667 ,  0.02103652,  0.01982814,  0.0183556 ,
                0.01663352,  0.01467847,  0.01250933,  0.01014766,  0.00761818,
                0.00494959,  0.00217557, -0.00066369, -0.00351907, -0.00632877,
                -0.00901202, -0.01145764, -0.01350184, -0.01487876, -0.01508707,
                -0.01287165,  0.        ])

Nx = x.shape[0]

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

# z-direction coordinates (span)
Nz = 51
z_max = 5
z = np.linspace(0, z_max, Nz)

# build 3D coordinate array
XYZ = np.zeros((Nx, Nz, 3))
XYZ[:, :, 0] = x[:, np.newaxis]
XYZ[:, :, 1] = y[:, np.newaxis]
XYZ[:, :, 2] = z[np.newaxis, :]


# apply spanwise-varying twist - max 45deg at tip
theta_twist = -np.pi/4

XYZ[:, :, 0] = (np.cos(theta_twist*XYZ[:, :, 2]/z_max)*XYZ[:, :, 0]
                - np.sin(theta_twist*XYZ[:, :, 2]/z_max)*XYZ[:, :, 1])
XYZ[:, :, 1] = (np.sin(theta_twist*XYZ[:, :, 2]/z_max)*XYZ[:, :, 0]
                + np.cos(theta_twist*XYZ[:, :, 2]/z_max)*XYZ[:, :, 1])


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# pass the mesh points to the StructuredGrid constructor, plot mesh points for
# visual inspection

airfoil = pv.StructuredGrid(XYZ[:, :, 0], XYZ[:, :, 1], XYZ[:, :, 2])

airfoil.plot(background='w', show_edges=True, pbr=True, metallic=1.0, roughness=0.4)


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Use pyvista to obtain surface normal vectors, plot normal vectors for
# inspection

# extract PolyData representation of mesh
airfoil2 = airfoil.extract_surface()

# extract normal vectors as numpy array
airfoil2 = airfoil2.compute_normals(flip_normals=True)
normals = (airfoil2.point_normals).reshape(XYZ.shape)

# plot mesh with normal vectors and save a screenshot
plotter = pv.Plotter()
plotter.add_mesh(airfoil2, show_edges=True)
plotter.add_arrows(airfoil2.points, airfoil2.point_normals, mag=0.02)

plotter.show_bounds(bounds=[0., 1, 0., 1, 0, 5], grid='front', location='outer',
                    all_edges=True)
plotter.show(auto_close=False)
# plotter.show(screenshot='my_normal_vectors.png')

# %%
# Build pywopwop patch object with structured grid geometry

# create empty instance of Geometry Patch File
rotorblade_geom = PWW.GeometryPatchFile()

# add header info
rotorblade_geom.geom_type = 'geometry'
rotorblade_geom.is_structured = 'y'
rotorblade_geom.time_type = 'constant'
rotorblade_geom.centered_type = 'node'
rotorblade_geom.float_type = 'single'
rotorblade_geom.iblank_included = 'n'

rotorblade_geom.units_string = 'Pa'
rotorblade_geom.comment_string = 'Test geometry file for single rotor blade'


# add rotor blade geometry and normals to patch obj
rotorblade_geom.addStructuredConstantZone('Rotor blade', XYZ.T, normals.T)

# print info
rotorblade_geom.print_info()

# create PSU-WOPWOP patch file
rotorblade_geom_filename = 'NACA4812_patch.dat'
rotorblade_geom.write_patch_file(rotorblade_geom_filename)
