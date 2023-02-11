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
import pathlib

from pywopwop._binary_readers_writers import read_block, write_block, read_int,\
    write_binary


# #############################################################################
# %% PSU-WOPWOP Sigma surface file converters - first approach
# #############################################################################

def process_sigma_geom_file(filename_geom, output_path='timesteps',
                            geometry_suffix='sigma_'):
    """
    Reads a multiple-timestep Sigma geometry (.x) file output from PSU-WOPWOP,
    and writes multiple single-timestep function (.x) files for opening in
    Paraview.

    Sigma single-timestep geometry files are enumerated sequentially as
    '[geometry_suffix][3-digit index].x' - e.g. 'sigma_000.x', 'sigma_001.x',
    etc.

    Parameters
    ----------
    filename_geom : str
        String containing the name of the geometry file - must include '.x'
        extension.

    output_path : str, optional
        Path where function will write the multiple single-timestep geometry
        files. Default is a new folder called 'timesteps'.

    geometry_suffix : str, optional
        Suffix to output filename. All output file names will consist of the
        given suffix plus a three-digit index. Default is 'sigma_'.

    Returns
    -------
    None


    Notes
    ----
    The geometry data is internally stored as a multidimensional list

        zones[nz][nx][nt][i, j]

    where:
    - nz: zone index
    - nx: coordinate index (0:x, 1:y, 2:z)
    - nt : time index
    """

    # *************************************************************************
    # check if output path exists; if it doesn't, create it. Also creates
    # parent paths, if necessary

    path = pathlib.Path(output_path)
    path.mkdir(parents=True, exist_ok=True)

    # *************************************************************************
    # read multiple-timestep geometry (.x) file

    with open(filename_geom, 'rb') as f:
        geom_data = f.read()

    # *************************************************************************
    # read file header

    # create lists of iMax, jMax, kMax, nVars for each block
    iMax_list = []
    jMax_list = []
    kMax_list = []          # 'k' is used as 'time' variable

    # read number of zones (i.e. independent meshes) in file
    Nzones = read_int(geom_data, 0)

    # read iMax, jMax, kMax for each zone
    for nz in range(Nzones):
        iMax_list.append(read_int(geom_data, 12*nz + 4))
        jMax_list.append(read_int(geom_data, 12*nz + 8))
        kMax_list.append(read_int(geom_data, 12*nz + 12))

    # end of header - set start index for beginning of geometry data
    start_index = 12*nz + 16

    # *************************************************************************
    # Read geometry data

    # create list of zones
    zones = []

    # for each zone...
    for nz in range(Nzones):

        # create list of timeseries for each coordinate
        timeseries_x = []
        timeseries_y = []
        timeseries_z = []

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # read and append 'x' coordinate data from each time step
        for it in range(kMax_list[nz]):
            block, start_index = read_block(geom_data, start_index, 1,
                                            iMax_list[nz], jMax_list[nz])
            timeseries_x.append(block)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # read and append 'y' coordinate data from each time step
        for it in range(kMax_list[nz]):
            block, start_index = read_block(geom_data, start_index, 1,
                                            iMax_list[nz], jMax_list[nz])
            timeseries_y.append(block)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # read and append 'z' coordinate data from each time step
        for it in range(kMax_list[nz]):
            block, start_index = read_block(geom_data, start_index, 1,
                                            iMax_list[nz], jMax_list[nz])
            timeseries_z.append(block)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        timeseries = [timeseries_x, timeseries_y, timeseries_z]

        # append current timestep list to block list
        zones.append(timeseries)


    # *************************************************************************
    # write multiple single-timestep geometry files

    # for each timestep...
    for nt in range(min(kMax_list)):

        with open(output_path + geometry_suffix + '{:03d}.x'.format(nt), 'wb') as file:

            # -----------------------------------------------------------------
            # write file header

            write_binary(file, Nzones)

            for nz in range(Nzones):
                write_binary(file, iMax_list[nz])
                write_binary(file, jMax_list[nz])
                write_binary(file, 1)

            # -----------------------------------------------------------------
            # write geometry data

            # for each block in current time step...
            for nz in range(Nzones):
                for nx in range(3):
                    write_block(file, zones[nz][nx][nt])
            # -----------------------------------------------------------------


# #############################################################################
def process_sigma_fn_file(filename_fn, filename_nam, output_path='timesteps',
                          function_suffix='sigma_'):
    """
    Reads a multiple-timestep Sigma function (.fn) file output from PSU-WOPWOP,
    and returns multiple single-timestep function (.fn) files for opening in
    Paraview.

    Sigma single-timestep function files are enumerated sequentially as
    '[function_suffix][3-digit index].fn' - e.g. 'sigma_000.fn',
    'sigma_001.fn', etc.


    Parameters
    ----------
    filename_fn : str
        Name of '.fn' file.

    filename_nam: str
        Name of '.nam' file.

    output_path: str, optional
        Path where function will write the multiple single-timestep function
        files. Default is a new folder called 'timesteps'.

    function_suffix : str, optional
        Suffix to output filename. All output file names will consist of the
        given suffix plus a three-digit index. Default is 'sigma_'.

    Returns
    ------
    sourcetime: (n_sourcetime,) array_like
        1D numpy array containing source times of 1st zone.


    Notes
    -----
    The return variable 'sourcetime' is used to tell Paraview what time (in
    seconds) corresponds to each time step. However, these values are taken
    from the first zone only, and thus might not correspond very accurately to
    the other zones.

    The function data is internally stored as a multidimensional list

        zones[nz][nvar][nt][i, j]

    where:
    - nz: zone index
    - nvar: variable index
    - nt : time index.

    """

    # *************************************************************************
    # check if output path exists; if it doesn't, create it. Also creates
    # parent paths, if necessary

    path = pathlib.Path(output_path)
    path.mkdir(parents=True, exist_ok=True)

    # *************************************************************************
    # read function (.fn) file
    with open(filename_fn, 'rb') as f:
        function_data = f.read()

    # *************************************************************************
    # read file header

    # read number of zones (i.e. independent meshes) in file
    Nzones = read_int(function_data, 0)

    # create lists of iMax, jMax, kMax, nVars for each block
    iMax_list = []
    jMax_list = []
    kMax_list = []          # 'k' is used as 'time' variable
    nVars_list = []

    for nz in range(Nzones):
        iMax_list.append(read_int(function_data, 16*nz + 4))
        jMax_list.append(read_int(function_data, 16*nz + 8))
        kMax_list.append(read_int(function_data, 16*nz + 12))
        nVars_list.append(read_int(function_data, 16*nz + 16))

    # end of header - set start index for beginning of function data
    start_index = 16*nz + 20

    # *************************************************************************
    # read function data

    # create list of zones
    zones = []

    # for each zone...
    for nz in range(Nzones):

        # create list of variables for current zone
        var_list = []

        # for each variable in current zone...
        for nvar in range(nVars_list[nz]):

            # create list of time steps for current variable
            timesteps = []

            # for each time step...
            for it in range(kMax_list[nz]):

                # read and append data from current time step
                #  ---> Data is always read as 1-dim; vector components are
                #       read as separate variables
                block, start_index = read_block(function_data, start_index,
                                                1, iMax_list[nz], jMax_list[nz])
                timesteps.append(block)

            # append current timestep list to variable list
            var_list.append(timesteps)

        # append current variable list to zone list
        zones.append(var_list)

    # ********************** Read source times *****************************
    # The source times of each zone are not identical; hence, determine what
    # source time ranges overlap with all other zones, and process only these

    # *************************************************************************
    # OLD CODE
    #
    # # extract vector of source times of first zone
    # n_sourcetime = min(kMax_list)
    #
    # sourcetime = np.zeros(n_sourcetime)
    # for it in range(n_sourcetime):
    #     sourcetime[it] = zones[0][0][it][0,0]
    #
    # *************************************************************************

    # create array of source times vs zones
    n_sourcetime = min(kMax_list)
    sourcetime = np.zeros((Nzones, n_sourcetime))
    for nz in range(Nzones):
        for nt in range(n_sourcetime):
            sourcetime[nz, nt] = zones[nz][0][nt][0,0]

    # calculate dt for source time (identical for all zones)
    dt = np.diff(sourcetime[0,:])[0]

    # find time shift per zone (vs zone 1), in number of 'source time' samples
    shift = np.round((sourcetime[:, 0] - sourcetime[0,0])/dt).astype(int)

    # calculate start and end indices for zones' source time
    start_nt = np.max(shift) - shift
    end_nt = n_sourcetime - shift

    # create vector of source times for zone 1 that overlap with all others
    sourcetime_final = sourcetime[0, start_nt[0]:end_nt[0]]

    # *************************************************************************
    # write multiple single-timestep function files

    # For each time step...
    for nt in range(sourcetime_final.shape[0]):

        # Create new "sigma_{:03d}.fn" files containing function data
        with open(output_path + function_suffix + '{:03d}.fn'.format(nt), 'wb') as file:

            # -----------------------------------------------------------------
            # write file header
            write_binary(file, Nzones)

            for nz in range(Nzones):
                write_binary(file, iMax_list[nz])
                write_binary(file, jMax_list[nz])
                write_binary(file, 1)
                write_binary(file, nVars_list[nz])

            # -----------------------------------------------------------------
            # write data

            # for each block in current time step...
            for nz in range(Nzones):

                # for each variable in current block...
                for nvar in range(nVars_list[nz]):
                    # write data within overlapping source time range
                    write_block(file, zones[nz][nvar][nt + start_nt[nz]])
            # -----------------------------------------------------------------

    return sourcetime_final


# #############################################################################
# %% PSU-WOPWOP Sigma surface file converters - second approach
# #############################################################################

def process_sigma_files(filename_geom, filename_fn, filename_nam,
                        output_path='timesteps/', geometry_suffix='sigma_',
                        function_suffix='sigma_',
                        p3d_filename='read_sigmasurfaces'):

    print("Begin processing Sigma files...")

    # extract names of Sigma variables from .nam file
    var_names = extract_sigma_var_names(filename_nam)
    print("Extracted sigma var names...")

    # read Sigma geometry file (multiple-timestep)
    zones_geom, geo_list = read_sigma_geom_file(filename_geom)

    # read Sigma function file (multiple-timestep)
    zones_fn, fn_list = read_sigma_fn_file(filename_fn)

    # process sigma source times
    sourcetime, start_nt = process_sigma_sourcetimes(zones_fn, fn_list)

    N_timesteps = sourcetime.shape[0]

    # write sigma geometry files (single-timestep)
    write_sigma_geom_files(zones_geom, geo_list, start_nt, N_timesteps,
                           output_path, geometry_suffix)

    # write Sigma function files (single-timestep)
    write_sigma_fn_files(zones_fn, fn_list, sourcetime, start_nt,
                         output_path, function_suffix)

    # write .p3d Paraview reader
    write_p3d_file(p3d_filename, output_path, sourcetime, var_names,
                   geometry_suffix, function_suffix)


# #############################################################################
def read_sigma_geom_file(filename_geom):
    """
    Reads a multiple-timestep Sigma geometry (.x) file, and returns a list of
    lists containing the geometry information.

    The geometry data returned is a multidimensional list of the form:

        zones[nz][nx][nt][i, j]

    where:
    - nz: zone index
    - nx: coordinate index (0:x, 1:y, 2:z)
    - nt : time index
    - i,j : structured mesh indices
    """

    # *************************************************************************
    # read multiple-timestep geometry (.x) file

    print('\nReading geometry file "' + filename_geom + '" ...')

    with open(filename_geom, 'rb') as f:
        geom_data = f.read()

    # *************************************************************************
    # Read file header

    # create lists of iMax, jMax, kMax, nVars for each block
    iMax_list = []
    jMax_list = []
    kMax_list = []          # 'k' is used as 'time' variable

    # read number of zones (i.e. independent meshes) in file
    Nzones = read_int(geom_data, 0)

    print('\tNzones = {}'.format(Nzones))

    # read iMax, jMax, kMax for each zone
    for nz in range(Nzones):
        iMax_list.append(read_int(geom_data, 12*nz + 4))
        jMax_list.append(read_int(geom_data, 12*nz + 8))
        kMax_list.append(read_int(geom_data, 12*nz + 12))

    # create output list with dims
    geo_list = [Nzones, iMax_list, jMax_list, kMax_list]

    # end of header - set start index for beginning of geometry data
    start_index = 12*nz + 16

    # *************************************************************************
    # Read geometry data

    # create list of zones
    zones_geom = []

    # for each zone...
    for nz in range(Nzones):
        print('\t\tReading zone {}...'.format(nz))

        # create list of timeseries for each coordinate
        timeseries_x = []
        timeseries_y = []
        timeseries_z = []

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # read and append 'x' coordinate data from each time step
        for it in range(kMax_list[nz]):
            block, start_index = read_block(geom_data, start_index, 1,
                                            iMax_list[nz], jMax_list[nz])
            timeseries_x.append(block)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # read and append 'y' coordinate data from each time step
        for it in range(kMax_list[nz]):
            block, start_index = read_block(geom_data, start_index, 1,
                                            iMax_list[nz], jMax_list[nz])
            timeseries_y.append(block)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # read and append 'z' coordinate data from each time step
        for it in range(kMax_list[nz]):
            block, start_index = read_block(geom_data, start_index, 1,
                                            iMax_list[nz], jMax_list[nz])
            timeseries_z.append(block)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        timeseries = [timeseries_x, timeseries_y, timeseries_z]

        # append current timestep list to block list
        zones_geom.append(timeseries)
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    print('\tFinished reading geometry file!')
    return zones_geom, geo_list


# #############################################################################
def write_sigma_geom_files(zones_geom, geo_list, start_nt, N_timesteps,
                           output_path='timesteps', geometry_suffix='sigma_'):
    """
    Receives a multidimensional list containing multiple-timestep Sigma surface
    geometry data, and writes multiple single-timestep function (.x) files for
    opening in Paraview.

    Sigma single-timestep geometry files are enumerated sequentially as
    '[geometry_suffix][3-digit index].x' - e.g. 'sigma_000.x', 'sigma_001.x',
    etc.

    Parameters
    ----------
    zones_geom : list
        Multidimensional list containing geometry data fom Sigma surface. See
        Notes below.

    geo_list : [Nzones, iMax_list, jMax_list, kMax_list] list
        TODO

    start_nt : int
        TODO

    N_timesteps : int
        TODO

    output_path : str, optional
        Path where function will write the multiple single-timestep geometry
        files. Default is a new folder called 'timesteps'.

    geometry_suffix : str, optional
        Suffix to output filename. All output file names will consist of the
        given suffix plus a three-digit index. Default is 'sigma_'.


    Returns
    -------
    None


    Notes
    ----
    The geometry data is input as a multidimensional list

        zones[nz][nx][nt][i, j]

    where:
    - nz: zone index
    - nx: coordinate index (0:x, 1:y, 2:z)
    - nt: time index
    """

    # *************************************************************************
    # check if output path exists; if it doesn't, create it. Also creates
    # parent paths, if necessary

    path = pathlib.Path(output_path)
    path.mkdir(parents=True, exist_ok=True)

    # *************************************************************************
    # parse input list for parameters

    print("\nWriting geometry files...")
    print("\tN_timesteps = {}".format(N_timesteps))

    Nzones, iMax_list, jMax_list, kMax_list = geo_list

    # *************************************************************************
    # Write multiple single-timestep geometry files

    # for each timestep...
    for nt in range(N_timesteps):
        print("\t\tNt = {}".format(nt))

        with open(output_path + geometry_suffix + '{:03d}.x'.format(nt), 'wb') as file:

            # -----------------------------------------------------------------
            # write file header

            write_binary(file, Nzones)

            for nz in range(Nzones):
                write_binary(file, iMax_list[nz])
                write_binary(file, jMax_list[nz])
                write_binary(file, 1)

            # -----------------------------------------------------------------
            # write geometry data

            # for each block in current time step...
            for nz in range(Nzones):
                for nx in range(3):
                    write_block(file, zones_geom[nz][nx][nt + start_nt[nz]])
            # -----------------------------------------------------------------

    print("\tFinished writing geometry files!")
    # *************************************************************************


# #############################################################################
def read_sigma_fn_file(filename_fn):

    # *************************************************************************
    # read function (.fn) file

    print('\nReading functional file "' + filename_fn + '" ...')

    with open(filename_fn, 'rb') as f:
        function_data = f.read()

    # *************************************************************************
    # read file header

    # read number of zones (i.e. independent meshes) in file
    Nzones = read_int(function_data, 0)

    print('\tNzones = {}'.format(Nzones))

    # create lists of iMax, jMax, kMax, nVars for each block
    iMax_list = []
    jMax_list = []
    kMax_list = []          # 'k' is used as 'time' variable
    nVars_list = []

    for nz in range(Nzones):
        iMax_list.append(read_int(function_data, 16*nz + 4))
        jMax_list.append(read_int(function_data, 16*nz + 8))
        kMax_list.append(read_int(function_data, 16*nz + 12))
        nVars_list.append(read_int(function_data, 16*nz + 16))


    fn_list = [Nzones, iMax_list, jMax_list, kMax_list, nVars_list]

    # end of header - set start index for beginning of function data
    start_index = 16*nz + 20

    # *************************************************************************
    # read function data

    # create list of zones
    zones_fn = []

    # for each zone...
    for nz in range(Nzones):

        print('\t\tReading zone {}...'.format(nz))

        # create list of variables for current zone
        var_list = []

        # for each variable in current zone...
        for nvar in range(nVars_list[nz]):

            # create list of time steps for current variable
            timesteps = []

            # for each time step...
            for it in range(kMax_list[nz]):

                # read and append data from current time step
                #  ---> Data is always read as 1-dim; vector components are
                #       read as separate variables
                block, start_index = read_block(function_data, start_index,
                                                1, iMax_list[nz], jMax_list[nz])
                timesteps.append(block)

            # append current timestep list to variable list
            var_list.append(timesteps)

        # append current variable list to zone list
        zones_fn.append(var_list)

    print("\tFinished reading functional file!")
    # *************************************************************************

    return zones_fn, fn_list


# #############################################################################
def process_sigma_sourcetimes(zones_fn, fn_list):

    # The source times of each zone are not identical; hence, determine what
    # source time ranges overlap with all other zones, and process only these

    Nzones, iMax_list, jMax_list, kMax_list, nVars_list = fn_list

    # create array of source times vs zones
    n_sourcetime = min(kMax_list)
    sourcetime = np.zeros((Nzones, n_sourcetime))
    for nz in range(Nzones):
        for nt in range(n_sourcetime):
            sourcetime[nz, nt] = zones_fn[nz][0][nt][0,0]

    # calculate dt for source time (identical for all zones)
    dt = np.diff(sourcetime[0,:])[0]

    # find time shift per zone (vs zone 1), in number of 'source time' samples
    shift = np.round((sourcetime[:, 0] - sourcetime[0,0])/dt).astype(int)

    # calculate start and end indices for zones' source time
    start_nt = np.max(shift) - shift
    end_nt = n_sourcetime - shift

    # create vector of source times for zone 1 that overlap with all others
    sourcetime_final = sourcetime[0, start_nt[0]:end_nt[0]]

    return sourcetime_final, start_nt


# #############################################################################
def write_sigma_fn_files(zones_fn, fn_list, sourcetime, start_nt,
                         output_path='timesteps', function_suffix='sigma_'):
    """
    fn_list : [Nzones, iMax_list, jMax_list, kMax_list, nVars_list]
    """

    print("\nWriting functional files...")
    # parse dimensions list
    Nzones, iMax_list, jMax_list, kMax_list, nVars_list = fn_list

    print("\tN_timesteps = {}".format(sourcetime.shape[0]))

    # *************************************************************************
    # write multiple single-timestep function files

    # For each time step...
    for nt in range(sourcetime.shape[0]):
        print("\t\tnt = {}".format(nt))

        # Create new "sigma_{:03d}.fn" files containing function data
        with open(output_path + function_suffix + '{:03d}.fn'.format(nt), 'wb') as file:

            # -----------------------------------------------------------------
            # write file header
            write_binary(file, Nzones)

            for nz in range(Nzones):
                write_binary(file, iMax_list[nz])
                write_binary(file, jMax_list[nz])
                write_binary(file, 1)
                write_binary(file, nVars_list[nz])

            # -----------------------------------------------------------------
            # write data

            # for each block in current time step...
            for nz in range(Nzones):

                # for each variable in current block...
                for nvar in range(nVars_list[nz]):
                    # write data within overlapping source time range
                    write_block(file, zones_fn[nz][nvar][nt + start_nt[nz]])
            # -----------------------------------------------------------------

    print("\tFinished writing functional files!")
    # *************************************************************************


# #############################################################################
# %% PSU-WOPWOP Sigma surface auxiliary functions
# #############################################################################



# def extract_sigma_var_names(filename_nam):
#    # *******************************************************************
#    # # --->>> MOVED TO _output_readers.read_namefile <<<---
#    # *******************************************************************
#     """
#     Parses a 'sigma.nam' file containing the variable names, and outputs a
#     list of these names.

#     Some vector components contain a semicolon in their name; if so, break
#     the name at the semicolon and keep just the 1st part.
#     """

#     var_names = []
#     with open(filename_nam, 'r') as file:
#         for line in file:
#             var_name = line.strip()

#             # check for semicolon
#             if ';' in var_name:
#                 var_name = var_name.split(';')[0]

#             var_names.append(var_name)

#     return var_names


# #############################################################################
def write_p3d_file(output_filename, output_path, sourcetime, var_names,
                   geometry_suffix='sigma_', function_suffix='sigma_'):
    """
    Writes .p3d reader file for Paraview.

    Parameters
    ----------
    output_filename : str
        Name of .p3d file (without .p3d extension).

    output_folder : str
        Path where function will write the .p3d file.

    sourcetime: (n_sourcetime,) array_like
        1D numpy array containing source times of 1st zone, obtained from
        'read_fn_file'.

    var_names: list of str
        List of sigma variable names, in order.

    geometry_suffix : str
        Suffix to single-timestep geometry filenames. Default is 'sigma_'.

    function_suffix : str, optional
        Suffix to single-timestep function filenames. Default is 'sigma_'.


    Returns
    -------
    None

    Notes
    -----
    Sigma files must be named as '[suffix][3-digit index].{x, fn}'. E.g.,
    'sigma_000.x' for geometry data of first timestep, 'sigma_001.x' for
    geometry data of 2nd timestep, etc.
    """

    p3d_header = ['{',
                  '\t"auto-detect-format" : true,',
                  '\t"filenames" : [']

    # *************************************************************************
    with open(output_path + output_filename + '.p3d', 'w') as file:

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # write p3d header
        for line in p3d_header:
            file.write(line+'\n')

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # write one line for each time step
        for nt in range(sourcetime.shape[0]):
            time_string = '\t\t{{"time" : {:.6f},'.format(sourcetime[nt])
            time_string += ' "xyz" : "' + geometry_suffix + '{:03d}.x",'.format(nt)
            time_string += ' "function" : "' + function_suffix + '{:03d}.fn"}},\n'.format(nt)
            file.write(time_string)

        file.write('\t\t],\n')

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # write sigma variables names
        function_string = '\t"function-names" : ['

        for i_name, function_name in enumerate(var_names):
            function_string += '"{}"'.format(function_name)

            if i_name < len(var_names)-1:
                function_string += ', '

        function_string += ']\n'
        file.write(function_string)
        file.write('}')
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    # *************************************************************************
