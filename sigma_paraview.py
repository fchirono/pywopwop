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

from readers_and_writers import read_block, write_block, read_int, write_binary
from consts_and_dicts import sigma_vars_dict


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


def read_fn_file(filename_list, output_path):
    """
    Reads a multiple-timestep Sigma function (.fn) file output from PSU-WOPWOP,
    and returns multiple single-timestep function (.fn) files for opening in
    Paraview.
    
    filename_list: list containing .x, .fn, .nam file names as strings.
    """
    
    # separate and parse filenames list
    for name in filename_list:
        
        file_extension = name.split('.')[1]
        
        # verify file name extension
        assert (file_extension in ['x', 'fn', 'nam']), \
            'File extension of {} in "filename_list" not recognized!'.format(name)
        
        if file_extension == 'x':
            filename_geom = name
        
        elif file_extension == 'fn':
            filename_function = name
    
        elif file_extension == 'nam':
            filename_names = name
            
            
    # extract variable names from .nam file
    var_names = extract_var_names(filename_names)
    
    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # read function (.fn) file
    
    with open(filename_function, 'rb') as f:
        function_data = f.read()
    
    
    # read number of blocks (i.e. independent meshes) in file
    Nblocks = read_int(function_data, 0)
    
    # create lists of iMax, jMax, kMax, nVars for each block
    iMax_list = []
    jMax_list = []
    kMax_list = []
    nVars_list = []
    
    # read header
    for n in range(Nblocks):
        
        iMax_list.append(read_int(function_data, 16*n + 4))
        jMax_list.append(read_int(function_data, 16*n + 8))
        kMax_list.append(read_int(function_data, 16*n + 12))     # 'time' var
        nVars_list.append(read_int(function_data, 16*n + 16))
    
    start_index = 16*n + 20

        
    blocks = []
    
    # for each block...
    for ib in range(Nblocks):
    
        # create list of variables for current block
        var_list = []
        
        # for each variable in current block...
        for ivar in range(nVars_list[ib]):
            
            # get number of dims of current variable (e.g. 1 for scalar data,
            # 3 for vector data)
            ndim_var = sigma_vars_dict[var_names[ivar]]
            
            # create list of time steps for current variable
            timesteps = []
            
            # for each time step...
            for it in range(kMax_list[ib]):
                
                # read and append data from current time step
                block, start_index = read_block(function_data, start_index,
                                                ndim_var, iMax_list[ib], jMax_list[ib])
                timesteps.append(block)
    
            # append current timestep list to variable list
            var_list.append(timesteps)
        
        # append current variable list to block list
        blocks.append(var_list)
    

    # access data as: blocks[ib][ivar][it][i, j]
    # - ib: block index
    # - ivar: variable index
    # - it : time index
    
    # extract vector of source times
    n_sourcetime = min(kMax_list)
    sourcetime = np.zeros(n_sourcetime)
    
    for it in range(n_sourcetime):
        sourcetime[it] = blocks[0][0][it][0,0]
    
    
    # Create new "sigma_{:03d}.fn" files containing function data per time step
    
    for nt in range(n_sourcetime):
        
        with open(output_path + 'sigma_{:03d}.fn'.format(nt), 'wb') as file:
            
            # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            # write file header
            write_binary(file, Nblocks)
            
            for ib in range(Nblocks):
                write_binary(file, iMax_list[ib])
                write_binary(file, jMax_list[ib])
                write_binary(file, 1)
                write_binary(file, nVars_list[ib])
            
            # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            # write data
            
            # for each block in current time step...
            for ib in range(Nblocks):
            
                # for each variable in current block...
                for ivar in range(nVars_list[ib]):
                    write_block(file, blocks[ib][ivar][nt])
            # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-     