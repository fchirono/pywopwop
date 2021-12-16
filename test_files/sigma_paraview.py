# -*- coding: utf-8 -*-
"""
pywopwop - https://github.com/fchirono/pywopwop


Parses a group of multiple-timesteps PSU-WOPWOP sigma surface files (.x, .nam,
.fn), and breaks them apart into separate geometry and function files per time
step for reading into Paraview.

For vector data, each component (x,y,z) is output as a separate variable. These
must be built into a vector array using, for example, the Python Calculator in
Paraview [1].

For example, for three velocity components 'Vx', 'Vy', 'Vz':
    - Load the .p3d file in Paraview;

    - Add a new filter via 'Filters' - 'Alphabetical' - 'Python Calculator';

    - In Python Calculator, set the 'Expression' to 'make_vector(Vx, Vy, Vz)'
    and the 'Array Name' to e.g. 'Velocity'. Click 'Apply'.

You should now have a new variable called 'Velocity' available for selection.
You can now choose to color the visualization by, e.g., vector magnitude in the
'Coloring' section.


To visualize the vectors as arrows:
    - Add a new glyph via 'Filters' - 'Alphabetical' - 'Glyph';

    - In 'Glyph Source' - 'Glyph Type', select 'Arrow';

    - In 'Orientation' - 'Orientation Array', select your vector - e.g. 'V'.


For a clearer visualization, you might need to manually adjust 'Scale' -
'Scale Factor'. Try playing with 'Scale' - 'Scale Array', and choose an
informative coloring scheme in 'Coloring' - e.g., vector magnitude.

For efficiency, you probably only need to visualize every few points. In
'Masking' - 'Glyph Mode', select 'Every Nth point' and in 'Stride' select
a number between, e.g., 2 and 20.


Author:
    Fabio Casagrande Hirono
    Dec 2021


References:
    [1] https://discourse.paraview.org/t/how-do-i-use-three-different-arrays-as-components-of-one-vector-field/4637
"""


import pywopwop as PWW


# %% #######################################################################
# define subfolder to analyse - run from within 'PSU-WOPWOP_v3.4.4' folder

# path where sigma files (.x, .nam, .fn) are located
path_to_sigma = './case1/'

# file names
filename_nam = path_to_sigma + '/sigma.nam'
filename_geom = path_to_sigma + '/sigma.x'
filename_fn = path_to_sigma + '/sigma.fn'

# path for file output - here we create a new folder called 'timesteps' to
# store the newly created single-timestep (.x, .fn) and Paraview reader (.p3d)
# files
path_output = path_to_sigma + '/timesteps2/'


# %% #######################################################################
# Approach 1: merge zones' data by time step index (i.e. all zones' first step
# are used to create the first single-timestep output file). This leads to
# mismatched source times, as zones have different values for their source
# times.
#
# This section has not been refactored into a neat, single function call!

# # extract names of Sigma variables from .nam file
# var_names = PWW.extract_sigma_var_names(filename_nam)

# # process geometry (.x) file
# PWW.process_sigma_geom_file(filename_geom, path_output)

# # process function files, and obtain a vector of source times for zone 1
# source_time = PWW.process_sigma_fn_file(filename_fn, filename_nam, path_output)

# # write Paraview .p3d reader file
# p3d_filename = 'read_sigmasurfaces_case1'
# PWW.write_p3d_file(p3d_filename, path_output, source_time, var_names)

# %% #######################################################################
# Approach 2: align zones to their nearest source time steps, discarding steps
# that only exist for some zones. This reduces the temporal/spatial
# misalignment and leads to nicer visualisations in Paraview.

PWW.process_sigma_files(filename_geom, filename_fn, filename_nam,
                        path_output)
