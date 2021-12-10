"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    --> PSU-WOPWOP Sigma surface file converter for visualization with
        Paraview.

Author:
    Fabio Casagrande Hirono
    Dec 2021
"""


import numpy as np

from readers_and_writers import initial_check, read_block, write_block, \
    read_IBLANKblock, read_int, read_float, write_binary, write_string, \
    read_string


def extract_var_names(nam_filename):
    """
    Parses a 'sigma.nam' file containing the variable names, and outputs a
    list of these names.
    """
    
    names = []
    with open(nam_filename, 'r') as file:
        for line in file:
            names.append(line.strip())
    
    return names


def read_fn_header(filename):
    """
    Reads the header of a Sigma function file (.fn) output from PSU-WOPWOP.
    """
    
    with open(filename, 'rb') as f:
        bytes_data = f.read()
    
    
    # read number of blocks (i.e. independent meshes) in file
    Nblocks = read_int(bytes_data, 0)
    
    # create lists of iMax, jMax, kMax, nVars for each block
    iMax_list = []
    jMax_list = []
    kMax_list = []
    nVars_list = []
    
    # read header
    for n in range(Nblocks):
        
        iMax_list.append(read_int(bytes_data, 16*n + 4))
        jMax_list.append(read_int(bytes_data, 16*n + 8))
        kMax_list.append(read_int(bytes_data, 16*n + 12))     # 'time' var
        nVars_list.append(read_int(bytes_data, 16*n + 16))
    
    start_index = 16*n + 20

    return (iMax_list, jMax_list, kMax_list, nVars_list, start_index)

