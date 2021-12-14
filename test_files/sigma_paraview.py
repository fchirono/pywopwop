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
    and the 'Array Name' to e.g. 'Velocity'. You can now choose to color the
    visualization by, e.g., vector magnitude in the 'Coloring' section.


To visualize the vectors as arrows:
    - Add a new glyph via 'Filters' - 'Alphabetical' - 'Glyph';
    
    - In 'Glyph Source' - 'Glyph Type', select 'Arrow';
    
    - In 'Orientation' - 'Orientation Array', select your vector - e.g. 'V';
    
    - In 'Scale' - 'Scale Array', select 'No scale array'. A manual 'Scale
    Factor' can be adjusted here for a better visualization;
    
    - For efficiency, you probably only need to visualize every few points. In
    'Masking' - 'Glyph Mode', select 'Every Nth point' and in 'Stride' select
    a number between, e.g., 2 and 20;
    
    - Pick an informative coloring scheme in 'Coloring' - e.g., vector magnitude.
    

Author:
    Fabio Casagrande Hirono
    Dec 2021


References:
    [1] https://discourse.paraview.org/t/how-do-i-use-three-different-arrays-as-components-of-one-vector-field/4637
"""


import pywopwop as PWW


# %% #######################################################################
# define subfolder to analyse - run from within 'PSU-WOPWOP_v3.4.4' folder

path_to_sigma = './case1/'

filename_nam = path_to_sigma + '/sigma.nam'
filename_geom = path_to_sigma + '/sigma.x'
filename_fn = path_to_sigma + '/sigma.fn'

# path for file output
path_output = path_to_sigma + '/timesteps/'


# %% #######################################################################

# extract names of Sigma variables from .nam file
var_names = PWW.extract_sigma_var_names(filename_nam)

# parse geometry (.x) file
PWW.process_sigma_geom_file(filename_geom, path_output)

# parse function files, and obtain a vector of source times
source_time = PWW.process_sigma_fn_file(filename_fn, filename_nam, path_output)

# write Paraview .p3d reader file
p3d_filename = 'read_sigmasurfaces_case1'
PWW.write_p3d_file(p3d_filename, path_output, source_time, var_names)

