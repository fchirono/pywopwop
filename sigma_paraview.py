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


def parse_filename_list(filename_list):
    # separate and parse list of filenames
    
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

    return filename_geom, filename_function, filename_names


def extract_var_names(nam_filename):
    """
    Parses a 'sigma.nam' file containing the variable names, and outputs a
    list of these names.
    """
    
    var_names = []
    with open(nam_filename, 'r') as file:
        for line in file:
            var_names.append(line.strip())
    
    return var_names


def read_fn_file(filename_function, filename_names, output_path):
    """
    Reads a multiple-timestep Sigma function (.fn) file output from PSU-WOPWOP,
    and returns multiple single-timestep function (.fn) files for opening in
    Paraview.
    
    filename_list: list containing .x, .fn, .nam file names as strings.
    """
    
            
            
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


def write_p3d_file(output_folder, output_filename, source_time, var_names):
    """
    Writes .p3d reader file for Paraview.
    
    Sigma files must be named sequentially: e.g. 'sigma_000.x' for geometry
    data of first timestep, 'sigma_001.x' for geometry data of 2nd timestep, etc.
    
    output_folder: where to save .p3d file
    output_filename: name of .p3d file (without .p3d extension)
    source_time: 1D numpy array of source times (taken from 1st block -
                                                 doesn't necessarily matches other blocks)
    var_names: list of sigma variable names
    """
    
    p3d_header = ['{',
                  '\t"auto-detect-format" : true,',
                  '\t"filenames" : [']
    
    with open(output_folder + output_filename + '.p3d', 'w') as file:
    
        # write p3d header
        for line in p3d_header:
            file.write(line+'\n')
    
        # write one line for each time step
        for ti in range(source_time.shape[0]):
            time_string = '\t\t{{"time" : {:.6f},'.format(source_time[ti])
            time_string += ' "xyz" : "sigma_{:03d}.x",'.format(ti)
            time_string += ' "function" : "sigma_{:03d}.fn"}},\n'.format(ti)
            file.write(time_string)
    
        file.write('\t\t],\n')
    
        # write sigma variables names
        function_string = '\t"function-names" : ['
        
        for i_name, function_name in enumerate(var_names):
            function_string += '"{}"'.format(function_name)
        
            if i_name < len(var_names)-1:
                function_string += ', '
        
        function_string += ']\n'
        file.write(function_string)
        file.write('}')
