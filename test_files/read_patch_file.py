# -*- coding: utf-8 -*-
"""
Reads PSU-WOPWOP binary patch (.dat) files using pywopwop

File contents can be compared in Cygwin using:
    cmp ./path-to/file1.dat ./path-to/file2.dat

Author:
    Fabio Casagrande Hirono
    Aug 2021
"""


import pywopwop as PWW
import matplotlib.pyplot as plt


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Opens patch files for inspection

# case 1
path = '../PSU-WOPWOP_v3.4.4/case1/'
filename = 'gyrodyne.dat'

# # case 2
# path = '../PSU-WOPWOP_v3.4.4/case2/'
# filename = 'reduced_geometry.dat'

# # case 5
# path = '../PSU-WOPWOP_v3.4.4/case5/'
# filename = 'constGeo_short.dat'


# create instance of Geometry Patch File and read file
case_patchfile = PWW.GeometryPatchFile(path+filename)

# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# writes content to different geometry patch file

filename2 = 'MyNewFile.dat'

case_patchfile.write_patch_file(filename2)


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Plot mesh points and normal vectors

from matplotlib import colors
from matplotlib import cm

import numpy as np

cmap = cm.get_cmap('jet')

# normalize pressure values
norm = colors.Normalize(vmin=0, vmax=4)

m_list =['o', '*', 's', 'p', '^']

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

for nz in range(case_patchfile.n_zones):
    # for nz in [0,]:

    # plot surface points
    XYZ = np.copy(case_patchfile.zones[nz].XYZ_coord)
    ax.scatter(XYZ[0], XYZ[1], XYZ[2], marker='o', color=cmap(norm(nz)))

    # ax.plot_surface(XYZ[0], XYZ[1], XYZ[2], linewidth=1.,
    #                 facecolors=cmap(norm(nz*np.ones((XYZ.shape[1], XYZ.shape[2])))))


    # plot normal vectors
    normal = np.copy(case_patchfile.zones[nz].normal_coord)
    # normal *= 1/np.linalg.norm(normal, 2, axis=0)     # normalize to unit magnitude
    ax.scatter(normal[0]+XYZ[0], normal[1]+XYZ[1], normal[2]+XYZ[2], marker='^',
                color=cmap(norm(nz)+0.1))


ax.set_title('XYZ Coords')


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()

