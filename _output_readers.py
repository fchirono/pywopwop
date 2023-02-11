"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    --> PSU-WOPWOP output file readers

Author:
    Fabio Casagrande Hirono
    Nov 2022
"""


import numpy as np


from pywopwop._consts_and_dicts import VALUE_LENGTH
from pywopwop._binary_readers_writers import read_block, read_int


# #############################################################################
# %% Readers for observer grid output files
# #############################################################################


def read_geometry_obs_grid(filename):
    """
    Reads PSU-WOPWOP output observer grid coordinate files.

    File is PLOT3D-like, contains a single zone, and the 'k' dimension holds
    the time information.

    Output is (x, y, z) grid with dimensions (3, Nt, iMax, jMax)
    """

    with open(filename, 'rb') as f:
        geom_data = f.read()

    # spherical grid: nbTheta
    iMax = read_int(geom_data, 0)

    # spherical grid: nbPsi
    jMax = read_int(geom_data, VALUE_LENGTH)

    # number of timesteps ('k' dimension)
    Nt = read_int(geom_data, 2*VALUE_LENGTH)

    # read blocks of obs coordinates
    start_index = 3*VALUE_LENGTH

    # x, y, z are (Nt, iMax, jMax) blocks
    x, next_index = read_block(geom_data, start_index, Nt, iMax, jMax)
    y, next_index = read_block(geom_data, next_index, Nt, iMax, jMax)
    z, next_index = read_block(geom_data, next_index, Nt, iMax, jMax)

    xyz = np.stack((x, y, z), axis=0)

    return xyz


def read_pressures_obs_grid(filename, remove_mean=True):
    """
    Reads PSU-WOPWOP output observer grid data files.

    File is PLOT3D-like, contains a single zone, and the 'k' dimension holds
    the time information.

    File can contain multiple variables, depending on how the PSU-WOPWOP case
    was set up. A typical example contains 4 variables: time, thickness noise,
    loading noise, and total noise. These are listed in the 'filename.nam'
    file.

    Output is a data array with dimensions (Nvar, Nt, iMax, jMax).

    If 'remove_mean' is true (default), each variable has its mean value
    removed before being output.
    """

    with open(filename, 'rb') as f:
        pressure_data = f.read()

    # spherical grid: nbTheta
    iMax = read_int(pressure_data, 0)

    # spherical grid: nbPsi
    jMax = read_int(pressure_data, VALUE_LENGTH)

    # number of timesteps ('k' dimension)
    Nt = read_int(pressure_data, 2*VALUE_LENGTH)

    # Number of variables to read
    #   --> e.g. time, thickness noise, loading noise, total noise
    Nvar = read_int(pressure_data, 3*VALUE_LENGTH)

    data = np.zeros((Nvar, Nt, iMax, jMax))

    next_index = 4*VALUE_LENGTH

    for n in range(Nvar):
        data[n, ...], next_index = read_block(pressure_data, next_index,
                                              Nt, iMax, jMax)

    # remove mean of each variable
    if remove_mean:
        for n in range(1, Nvar):
            mean = np.mean(data[n, ...], axis=0)
            data[n, ...] += -mean

    return data


# #############################################################################
# %% Readers for single-observer output files
# #############################################################################

def read_pressures_single_obs(filename, print_varnames=False):
    """
    Reads a PSU-WOPWOP output file ('.tec') containing the acoustic pressure(s)
    time history at a single observer, and returns the output as a dictionary.

    The file 'filename.tec' is in ASCII Tecplot 'POINT' format, and contains a
    list of variable names. These variable names are used as the keys for the
    output dictionary, while the values are one-dimensional Numpy arrays.

    A typical PSU-WOPWOP run might contains the following variables:
        - 'Observer Time (s)'
        - 'Thickness Acoustic Pressure (Pa)'
        - 'Loading Acoustic Pressure (Pa)'
        - 'Total Acoustic Pressure (Pa)'


    Parameters
    ----------
    filename : str
        String containing the file to be read.

    print_varnames : bool, optional
        Boolean flag defining whether to print the variables' names found in
        the file. Default is False.

    Returns
    -------
    output : dict
        Dictionary containing the variables

    """

    # First four rows do not contain numerical data
    rows_to_skip = 4

    with open(filename, 'r') as f:
        lines = f.readlines()

    assert "PSU-WOPWOP" in lines[0], \
        "Error when reading file {} header: this might not be a PSU-WOPWOP single-observer output file!"

    var_list = _retrieve_var_list(lines)
    values = np.loadtxt(filename, skiprows=rows_to_skip)

    output = {}
    for i, var in enumerate(var_list):
        output[var] = np.copy(values[:, i])

    if print_varnames:
        print('List of variables in {}:'.format(filename))
        for var in output.keys():
            print('\t{}'.format(var))

    return output


# #############################################################################
def _retrieve_var_list(lines):
    """
    Reads the 2nd line of a PSU-WOPWOP single-observer output file and outputs
    a list of 'clean' variable names.
    """

    import re

    # parse 2nd line for list of variables (ignore first few chars, split at
    # commas)
    raw_varlist = lines[1][10:].split(',')

    varlist = []

    for v in raw_varlist:
        # use regular expressions to obtain substring between double quotes
        varlist.append(re.findall('"([^"]*)"', v)[0])

    return varlist


# #############################################################################
# %% Auxiliary functions
# #############################################################################

def read_namefile(filename):
    """
    Parses a name file ('.nam') containing the variable names, and outputs a
    list of these names.

    When reading Sigma Surfaces, some vector components contain a semicolon in
    their name; if so, break the name at the semicolon and keep just the 1st
    part.
    """

    var_names = []
    with open(filename, 'r') as file:
        for line in file:
            var_name = line.strip()

            # check for semicolon
            if ';' in var_name:
                var_name = var_name.split(';')[0]

            var_names.append(var_name)

    return var_names


# #############################################################################
def interp_fs(t_original, x_original, fs):
    """
    Interpolates an array of non-uniformly sampled values 'x_original', with
    samples acquired at time instants 't_original', over a uniformly sampled
    array 't_new' sampled at 'fs' Hz.

    Parameters
    ----------
    t_original : (N_samples,) array_like
        Array of non-uniform sampling times.

    x_original : (N_samples,) array_like
        Array containing non-uniformly sampled values.

    fs : int
        Sampling frequency (in Hz) for the new uniformly sampled time vector.

    Returns
    -------
    t_unif : (N_new,) array_like
        Array of uniformly sampled time instants.

    x_unif : (N_new,) array_like
        Array of uniformly sampled, interpolated values.

    Notes
    -----
    The new time samples are all contained inside the range of original time
    samples 't_original'.
    """

    from scipy.interpolate import interp1d

    # find min, max time values within original time samples
    nt_min = np.ceil(t_original[0]*fs).astype(int)
    nt_max = np.floor(t_original[-1]*fs).astype(int)

    t_min = nt_min/fs
    t_max = nt_max/fs

    # number of time samples in uniform array
    Nt = ((t_max-t_min)*fs).astype(int)

    t_unif = np.linspace(t_min, t_max, Nt)

    x_interpolator = interp1d(t_original, x_original)

    return t_unif, x_interpolator(t_unif)
