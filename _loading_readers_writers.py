# -*- coding: utf-8 -*-
"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    This file contains methods for reading and writing PSU-WOPWOP loading
    files. These methods belong to the 'PWWPatch' class, defined in
    '__init__.py'.

Author:
    Fabio Casagrande Hirono
    Nov 2022
"""

import numpy as np

from ._consts_and_dicts import MAGICNUMBER, ENDIANNESS, VALUE_LENGTH, \
    reverse_dict, structured_dict, loading_time_dict, structured_header_length, \
    centered_dict, loading_data_dict, ref_frame_dict, float_dict

from ._binary_readers_writers import initial_check, read_block, write_block, \
    read_int, read_float, write_binary, write_string, read_string


# ##########################################################################
# %%PSU-WOPWOP loading file readers
# ##########################################################################

def _read_loading_header(self, loading_filename):
    """
    Reads PSU-WOPWOP loading functional header info from 'loading_filename'.
    """

    # do initial check for 'magic number' and endianness
    file_endianness = initial_check(loading_filename)
    assert file_endianness == ENDIANNESS, \
        "File endianness does not match pywopwop default!"

    # open file and read binary content
    with open(loading_filename, 'rb') as f:
        bytes_data = f.read()

    # check version number - must be '1'and '0'
    self.version_number_major = read_int(bytes_data, 4)
    self.version_number_minor = read_int(bytes_data, 8)

    assert ((self.version_number_major == 1) and (self.version_number_minor == 0)), \
        'File version is not v1.0!'

    # read comments string (1024 bytes, starting at index 12)
    loading_comment = read_string(bytes_data, 12, 1024)
    self.set_loading_comment(loading_comment)


    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # read format string (10 ints, 40 bytes, starting at index 1036),
    # assert loading file info vs. previously read from geometry file

    self.loading_format_string = []
    for n in range(10):
        self.loading_format_string.append(read_int(bytes_data, 1036 + 4*n))

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Verify file type description

    # --> loading_format_string[0] is '2' to indicate functional data file
    assert (self.loading_format_string[0] == 2), \
        "Format string does not start with '2', not a PSU-WOPWOP functional file!"

    assert (self.n_zones == self.loading_format_string[1]), \
        "Number of zones in loading file does not match geometry file!"

    assert (self.is_structured
            == reverse_dict(structured_dict, self.loading_format_string[2])), \
        "Loading file 'is_structured' property does not match!"

    self.loading_time_type = \
        reverse_dict(loading_time_dict, self.loading_format_string[3])

    assert (self.centered_type
            == reverse_dict(centered_dict, self.loading_format_string[4])), \
        "Loading file 'centered_type' property does not match!"

    self.loading_data_type  = \
        reverse_dict(loading_data_dict, self.loading_format_string[5])

    self.loading_ref_frame  = \
        reverse_dict(ref_frame_dict, self.loading_format_string[6])

    assert (self.float_type
            == reverse_dict(float_dict, self.loading_format_string[7])), \
        "Loading file 'float_type' property does not match!"

    # --> loading_format_string[8] is reserved for future use, and must be '0' in this version
    # --> loading_format_string[9] is reserved for future use, and must be '0' in this version


    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # zone specification

    # read number of zones containing functional data (4 bytes)
    self.n_zones_with_loading_data = read_int(bytes_data, 1076)

    # create list of zones containing functional data
    # -->>  negative numbers indicate zones for which WOPWOP should *NOT*
    #       calculate thickness noise - e.g. loading patch
    #
    # -->> PSU_WOPWOP zone list indices are one-based - subtract one to
    #       obtain zero-based (Python) indices
    self.zones_with_loading_data = []

    for z in range(self.n_zones_with_loading_data):
        self.zones_with_loading_data.append(read_int(bytes_data, 1080 + z*4) - 1)


    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # read header info

    zone_info_start = 1080 + self.n_zones_with_loading_data*4

    # version 1.0
    if self.version_number_minor == 0:

        # --------------------------------------------------------------
        # structured loading
        if self.is_structured == True:

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if self.loading_time_type == 'constant':

                for i, nz in enumerate(self.zones_with_loading_data):

                    # remove sign from zone index
                    nz = abs(nz)

                    # get handle to existing zone in list
                    zone = self.zones[nz]

                    # set loading header length
                    zone.loading_header_length = \
                        structured_header_length[self.loading_time_type]

                    # read loading zone name
                    name = read_string(bytes_data,
                                       zone_info_start + i*zone.loading_header_length, 32)
                    zone._set_string(name, 'loading_name', 32)

                    # assert iMax and jMax match
                    iMax_fromfile = read_int(bytes_data,
                                             zone_info_start + 32 + i*zone.loading_header_length)
                    jMax_fromfile = read_int(bytes_data,
                                             zone_info_start + 36 + i*zone.loading_header_length)
                    assert ((zone.iMax == iMax_fromfile) and (zone.jMax == jMax_fromfile)), \
                        "(iMax, jMax) in Zone {} from loading file don't match existing values in PWWPatch instance!".format(nz)

                    # set loading data flag
                    zone.has_loading_data = True

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            elif self.loading_time_type == 'aperiodic':

                for i, nz in enumerate(self.zones_with_loading_data):

                    # remove sign from zone index
                    nz = abs(nz)

                    # get handle to existing zone in list
                    zone = self.zones[nz]

                    # set loading header length
                    zone.loading_header_length = \
                        structured_header_length[self.loading_time_type]

                    # read loading zone name
                    name = read_string(bytes_data,
                                       zone_info_start + i*zone.loading_header_length, 32)
                    zone._set_string(name, 'loading_name', 32)

                    # .........................................................
                    # read number of timesteps
                    loading_Nt = read_int(bytes_data,
                                          zone_info_start + 32 + i*zone.loading_header_length)

                    # check if 'Nt' attribute already exists (e.g. from aperiodic geometry)
                    if hasattr(self, 'Nt'):
                        # check for match
                        assert loading_Nt == self.Nt, \
                            "Number of timesteps in Zone {} loading data does not match existing PWWPatch instance 'Nt'!".format(nz)
                    else:
                        # store 'Nt' in PWWPatch
                        self.Nt = loading_Nt

                    # .........................................................
                    # assert iMax and jMax match
                    iMax_fromfile = read_int(bytes_data,
                                             zone_info_start + 36 + i*zone.loading_header_length)
                    jMax_fromfile = read_int(bytes_data,
                                             zone_info_start + 40 + i*zone.loading_header_length)

                    assert ((zone.iMax == iMax_fromfile) and (zone.jMax == jMax_fromfile)), \
                        "(iMax, jMax) in Zone {} from loading file don't match existing values in PWWPatch instance!".format(nz)

                    # set loading data flag
                    zone.has_loading_data = True

                # check if 'time_steps' attribute already exists (e.g. from aperiodic geometry)
                if not hasattr(self, 'time_steps'):
                    # store 'Nt' in PWWPatch
                    self.time_steps = np.zeros(self.Nt, dtype=np.float32)

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            else:
                # TODO: implement non-constant non-aperiodic loading
                raise NotImplementedError("Can't read loading data that is not constant nor aperiodic - not implemented yet!")
            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

        else:
            # TODO: implement non-structured loading
            raise NotImplementedError("Can't read non-structured loading data - not implemented yet!")

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # TODO: implement header reader for version 1.1
    elif self.version_number2 == '1':
        raise NotImplementedError("Can't read loading data file v1.1 header - not implemented yet!")


# ##########################################################################
def _read_loading_data(self, loading_filename):
    """
    Reads PSU-WOPWOP loading data from 'loading_filename'.
    """

    # open file and read binary content
    with open(loading_filename, 'rb') as f:
        bytes_data = f.read()

    # *************************************************************************
    # structured loading
    if self.is_structured == True:

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        if self.loading_time_type == 'constant':

            # start index for reading functional data
            # end of format string + zone specification + header
            field_start = (1076
                           + (1 + self.n_zones_with_loading_data)*4
                           + self.n_zones_with_loading_data*self.zones[0].loading_header_length)

            # ----------------------------------------------------------------
            # if data is surface pressure
            if self.loading_data_type == 'surf_pressure':

                # for each zone with loading data
                for nz in self.zones_with_loading_data:

                    # create empty numpy arrays for pressure data
                    pressures = np.zeros((self.zones[nz].iMax, self.zones[nz].jMax),
                                         dtype=np.float32)

                    self.zones[nz].add_StructuredConstantLoading(pressures, 'surf_pressure')

                    # read pressure data and next index
                    self.zones[nz].loading.pressures, field_start = \
                        read_block(bytes_data, field_start, 1,
                                   self.zones[nz].iMax,
                                   self.zones[nz].jMax)

            # ----------------------------------------------------------------
            # if data is loading vectors
            elif self.loading_data_type == 'surf_loading_vec':

                # for each zone containing data:
                for nz in self.zones_with_loading_data:

                    # create empty numpy arrays for pressure data
                    loading_vectors = np.zeros((3, self.zones[nz].iMax,
                                                self.zones[nz].jMax),
                                               dtype=np.float32)

                    self.zones[nz].add_StructuredConstantLoading(loading_vectors, 'surf_loading_vec')

                    # read pressure data and next index
                    self.zones[nz].loading.loading_vectors, field_start = \
                        read_block(bytes_data, field_start, 3,
                                   self.zones[nz].iMax,
                                   self.zones[nz].jMax)

            # ----------------------------------------------------------------
            elif self.loading_data_type == 'flow_params':

                # for each zone containing data:
                for nz in self.zones_with_loading_data:

                    # create empty numpy arrays for flow params
                    # (rho, rho*u, rho*v, rho*w, p')
                    flow_params = np.zeros((5, self.zones[nz].iMax,
                                            self.zones[nz].jMax),
                                           dtype=np.float32)

                    self.zones[nz].add_StructuredConstantLoading(flow_params, 'flow_params')

                    # read pressure data and next index
                    self.zones[nz].loading.flow_params, field_start = \
                        read_block(bytes_data, field_start, 5,
                                   self.zones[nz].iMax,
                                   self.zones[nz].jMax)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        elif self.loading_time_type == 'aperiodic':

            # start index for reading functional data
            # end of format string + zone specification + header
            field_start = (1076
                           + (1 + self.n_zones_with_loading_data)*VALUE_LENGTH
                           + self.n_zones_with_loading_data*self.zones[0].loading_header_length)

            # -----------------------------------------------------------------
            # if data is surface pressure
            if self.loading_data_type == 'surf_pressure':

                # for each zone with loading data, add empty numpy arrays for
                # pressure data and time steps
                for nz in self.zones_with_loading_data:

                    pressures = np.zeros((self.Nt, self.zones[nz].iMax,
                                          self.zones[nz].jMax),
                                         dtype=np.float32)

                    self.zones[nz].add_StructuredAperiodicLoading(pressures, 'surf_pressure')

                    # initialize 'time_steps' attribute in current Zone if it
                    # doesn't already exist (e.g. from reading aperiodic
                    # geometry);
                    if not hasattr(self.zones[nz], 'time_steps'):
                        self.zones[nz].time_steps = np.zeros((self.Nt,), dtype=np.float32)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # for each timestep...
                for nt in range(self.Nt):
                    # for each zone...
                    for nz in self.zones_with_loading_data:

                        # read current time value and next index
                        self.zones[nz].time_steps[nt], field_start = read_float(bytes_data, field_start)

                        # read pressure data and next index
                        self.zones[nz].loading.pressures[nt, :, :], field_start = \
                            read_block(bytes_data, field_start, 1,
                                       self.zones[nz].iMax,
                                       self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # TODO: assert all zones' time steps are identical!
                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

            # -----------------------------------------------------------------
            # if data is loading vectors
            elif self.loading_data_type == 'surf_loading_vec':

                # for each zone with loading data, add empty numpy arrays for
                # loading vector data and time steps
                for nz in self.zones_with_loading_data:

                    loading_vectors = np.zeros((self.Nt, 3,
                                                self.zones[nz].iMax,
                                                self.zones[nz].jMax),
                                               dtype=np.float32)

                    self.zones[nz].add_StructuredAperiodicLoading(loading_vectors,
                                                                  'surf_loading_vec')

                    # initialize 'time_steps' attribute in current Zone if it
                    # doesn't already exist (e.g. from reading aperiodic
                    # geometry);
                    if not hasattr(self.zones[nz], 'time_steps'):
                        self.zones[nz].time_steps = np.zeros((self.Nt,), dtype=np.float32)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # for each timestep...
                for nt in range(self.Nt):
                    # for each zone...
                    for nz in self.zones_with_loading_data:

                        # read current time value and next index
                        self.zones[nz].time_steps[nt], field_start = read_float(bytes_data, field_start)

                        # read loading vector data and next index
                        self.zones[nz].loading.loading_vectors[nt, :, :, :], field_start = \
                            read_block(bytes_data, field_start, 3,
                                       self.zones[nz].iMax,
                                       self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # TODO: assert all zones' time steps are identical!
                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

            # -----------------------------------------------------------------
            elif self.loading_data_type == 'flow_params':

                # for each zone with loading data, add empty numpy arrays for
                # flow parameter data and time steps
                for nz in self.zones_with_loading_data:

                    flow_params = np.zeros((self.Nt, 5,
                                            self.zones[nz].iMax,
                                            self.zones[nz].jMax),
                                           dtype=np.float32)

                    self.zones[nz].add_StructuredAperiodicLoading(flow_params,
                                                                  'flow_params')

                    # initialize 'time_steps' attribute in current Zone if it
                    # doesn't already exist (e.g. from reading aperiodic
                    # geometry);
                    if not hasattr(self.zones[nz], 'time_steps'):
                        self.zones[nz].time_steps = np.zeros((self.Nt,), dtype=np.float32)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # for each timestep...
                for nt in range(self.Nt):
                    # for each zone...
                    for nz in self.zones_with_loading_data:

                        # read current time value and next index
                        self.zones[nz].time_steps[nt], field_start = read_float(bytes_data, field_start)

                        # read flow parameter data and next index
                        self.zones[nz].loading.flow_params[nt, :, :, :], field_start = \
                            read_block(bytes_data, field_start, 5,
                                       self.zones[nz].iMax,
                                       self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # TODO: assert all zones' time steps are identical!
                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

            # -----------------------------------------------------------------
            # compare time_steps across all zones, ensure they are identical
            self.time_steps = np.copy(self.zones[0].time_steps)
            for nz in range(self.n_zones):
                assert np.allclose(self.zones[nz].time_steps, self.time_steps),\
                    "Error reading file {}: Zone {} time steps do not match Zone 0 time steps!".format(loading_filename, nz)
            # -----------------------------------------------------------------

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        else:
            # TODO: read non-constant, non-aperiodic loading data
            raise NotImplementedError("Can't read non-constant or non-aperiodic loading data - not implemented yet!")

    # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    else:
        # TODO: read non-structured loading data
        raise NotImplementedError("Can't read non-structured loading data - not implemented yet!")



# ##########################################################################
# %%PSU-WOPWOP loading file writers
# ##########################################################################

def _write_loading_header(self, loading_filename):
    """
    Writes PSU-WOPWOP functional file header to 'loading_filename'
    """

    # build format string with most recent format values
    self._build_loading_format_string()

    with open(loading_filename, 'wb') as file:

        # write 'magic number' 42 to first 4 bytes
        write_binary(file, MAGICNUMBER)

        # write version number
        write_binary(file, self.version_number_major)
        write_binary(file, self.version_number_minor)

        # write comments string (1024 bytes)
        write_string(file, self.loading_comment, 1024)

        # write format string (10 ints, 40 bytes)
        for n in range(10):
            write_binary(file, self.loading_format_string[n])

        # --------------------------------------------------------------
        # write zone specification

        # write number of zones with data
        write_binary(file, self.n_zones_with_loading_data)

        # write list of those zones - add one to create one-based
        # (PSU-WOPWOP) indices
        for z in self.zones_with_loading_data:
            write_binary(file, (z + 1))

        # --------------------------------------------------------------
        # write zone info

        if self.is_structured == True:

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if self.loading_time_type == 'constant':

                # for each zone containing data...
                for nz in self.zones_with_loading_data:

                    # write name (32-byte string)
                    write_string(file, self.zones[nz].loading_name, 32)

                    # write iMax and jMax (4 byte ints)
                    write_binary(file, self.zones[nz].iMax)
                    write_binary(file, self.zones[nz].jMax)

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            elif self.loading_time_type == 'aperiodic':

                # for each zone containing data...
                for nz in self.zones_with_loading_data:

                    # write name (32-byte string)
                    write_string(file, self.zones[nz].loading_name, 32)

                    # write number of timesteps
                    write_binary(file, self.Nt)

                    # write iMax and jMax (4 byte ints)
                    write_binary(file, self.zones[nz].iMax)
                    write_binary(file, self.zones[nz].jMax)
            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            else:
                # TODO: implement non-constant, non-aperiodic functional data header
                raise NotImplementedError("Can't write non-constant, non-aperiodic loading data header - not implemented yet!")

        # ------------------------------------------------------------
        else:
            # TODO: implement non-structured headers
            raise NotImplementedError("Can't write non-structured loading data header - not implemented yet!")

        # ------------------------------------------------------------



# ##########################################################################
def _write_loading_data(self, loading_filename):
    """
    Writes PSU-WOPWOP functional data to 'loading_filename'
    Must be immediately preceded by a call to
    '_write_loading_header(loading_filename)'!
    """

    # open file in append mode - no need to adjust index
    with open(loading_filename, 'ab') as file:

        # -----------------------------------------------------------
        # structured loading
        if self.is_structured == True:

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            # constant geometry
            if self.loading_time_type == 'constant':

                # ......................................................
                if self.loading_data_type == 'surf_pressure':

                    # for each zone...
                    for nz in self.zones_with_loading_data:
                        nz = abs(nz)
                        write_block(file, self.zones[nz].loading.pressures)

                # ......................................................
                elif self.loading_data_type == 'surf_loading_vec':

                    # for each zone...
                    for nz in self.zones_with_loading_data:
                        nz = abs(nz)
                        write_block(file, self.zones[nz].loading.loading_vectors)

                # ......................................................
                elif self.loading_data_type == 'flow_params':

                    # for each zone...
                    for nz in self.zones_with_loading_data:
                        nz = abs(nz)
                        write_block(file, self.zones[nz].loading.flow_params)

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            elif self.loading_time_type == 'aperiodic':

                # ......................................................
                if self.loading_data_type == 'surf_pressure':

                    # for each time step...
                    for nt, time in enumerate(self.time_steps):
                        # for each zone...
                        for nz in self.zones_with_loading_data:
                            nz = abs(nz)
                            write_binary(file, time)
                            write_block(file, self.zones[nz].loading.pressures[nt, :, :])

                # ......................................................
                elif self.loading_data_type == 'surf_loading_vec':

                    # for each time step...
                    for nt, time in enumerate(self.time_steps):
                        # for each zone...
                        for nz in self.zones_with_loading_data:
                            nz = abs(nz)
                            write_binary(file, time)
                            write_block(file, self.zones[nz].loading.loading_vectors[nt, :, :, :])

                # ......................................................
                elif self.loading_data_type == 'flow_params':

                    # for each time step...
                    for nt, time in enumerate(self.time_steps):
                        # for each zone...
                        for nz in self.zones_with_loading_data:
                            nz = abs(nz)
                            # write time value and block of flow params
                            write_binary(file, time)
                            write_block(file, self.zones[nz].loading.flow_params[nt, :, :, :])

            # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            else:
                # TODO: write non-constant, non-aperiodic loading data
                raise NotImplementedError("Can't write non-constant, non-aperiodic loading data - not implemented yet!")

        # -----------------------------------------------------------
        else:
            # TODO: write non-structured loading data
            raise NotImplementedError("Can't write non-structured loading data - not implemented yet!")


