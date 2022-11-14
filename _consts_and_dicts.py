"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    --> PSU-WOPWOP constants and dictionaries

Author:
    Fabio Casagrande Hirono
    Nov 2022
"""


# ##########################################################################
# %% PSU-WOPWOP program constants
# ##########################################################################

MAGICNUMBER = 42

# values are little endian, 4-byte, signed by default
ENDIANNESS = 'little'
VALUE_LENGTH = 4
IS_SIGNED = True

# Digit reserved for future use - must be '0' in this version
RESERVED_DIGIT = 0


# ##########################################################################
# %% PSU-WOPWOP program dictionaries and reverse_dict function
# ##########################################################################

def reverse_dict(my_dict, my_value):
    """
    Tries to find which key corresponds to a given value in a dictionary - i.e.
    reverses the dict lookup.

    This code assumes no two keys share the same value - as is the case with
    PSU-WOPWOP dicts used here.
    """

    return {value: key for key, value in my_dict.items()}[my_value]


# whether file is geometry file or node subset geometry file
geom_dict = {'geometry'     : 1,
             'node_subset'  :-1}

# whether file is structured or not
structured_dict = {True : 1,
                   False : 2}

# type of geometry (temporal behaviour)
# --> for loading files:
loading_time_dict = {'constant'                 :1,
                     'periodic'                 :2,
                     'aperiodic'                :3,
                     'mult_time_aperiodic'      :4}

# for geometry files: all of functional types, plus quasiperiodic types
geometry_time_dict = {**loading_time_dict,
                      'quasiperiodic'            :5,
                      'mult_time_quasiperiodc'   :6}

# whether normal vectors and areas are node-centered or face-centered
centered_dict = {'node' : 1,
                 'face' : 2}

# what type of loading data is contained in functional file
loading_data_dict = {'surf_pressure'    : 1,
                     'surf_loading_vec' : 2,
                     'flow_params'      : 3}

# reference frame for functional file (loading vectors and fluid momentum)
ref_frame_dict = {'ground_fixed': 1,
                  'mixed_frame' : 2,
                  'blade_fixed' : 3}

# whether floats are single or double precision
# --->>>    DOUBLES ARE NOT YET SUPPORTED! THIS PARAMETER IS
#           INTERNALY ASSUMED TO BE 1           <<<---
float_dict = {'single' : 1,
              'double' : 2}

# whether 'iblank' values are included with the geometry grid or not
iblank_dict = {True : 1,
               False: 0}

