# -*- coding: utf-8 -*-
"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    This file contains methods for reading and writing PSU-WOPWOP geometry
    files. These methods belong to the 'PWWPatch' class, defined in
    '__init__.py'.

Author:
    Fabio Casagrande Hirono
    Nov 2022
"""

import numpy as np

from ._consts_and_dicts import MAGICNUMBER, ENDIANNESS, reverse_dict, \
    geom_dict, structured_dict, geometry_time_dict, structured_header_length, \
    centered_dict, float_dict, iblank_dict

from ._binary_readers_writers import initial_check, read_block, write_block, \
    read_IBLANKblock, read_int, read_float, write_binary, write_string, \
    read_string

from ._zones import StructuredZone



# #############################################################################
# %%PSU-WOPWOP geometry file readers
# #############################################################################

def _read_geometry_header(self, geometry_filename):
    """
    Reads PSU-WOPWOP geometry header info from 'geometry_filename'.
    """

    # *************************************************************************
    # do initial check for 'magic number' and endianness
    file_endianness = initial_check(geometry_filename)
    assert file_endianness == ENDIANNESS, \
        "File endianness does not match pywopwop default!"

    # open fileand read binary content
    with open(geometry_filename, 'rb') as f:
        bytes_data = f.read()


    # check version number - must be '1'and '0'
    self.version_number_major = read_int(bytes_data, 4)
    self.version_number_minor = read_int(bytes_data, 8)

    assert ((self.version_number_major == 1)
            and (self.version_number_minor == 0)), \
        'File version is not v1.0!'

    # read units string (32 chars, starting at index 12)
    units_string = read_string(bytes_data, 12, 32)
    self.set_units_string(units_string)

    # read comments string (1024 bytes, starting at index 44)
    geom_comment = read_string(bytes_data, 44, 1024)
    self.set_geometry_comment(geom_comment)

    # *************************************************************************
    # read format string (8 ints, 32 bytes, starting at index 1068)
    self.geometry_format_string = []
    for n in range(8):
        self.geometry_format_string.append(read_int(bytes_data, 1068 + 4*n))

    # Populate file type description
    self.geometry_type      = reverse_dict(geom_dict,
                                           self.geometry_format_string[0])

    self.n_zones            = self.geometry_format_string[1]

    self.is_structured      = reverse_dict(structured_dict,
                                           self.geometry_format_string[2])

    self.geometry_time_type = reverse_dict(geometry_time_dict,
                                           self.geometry_format_string[3])

    self.centered_type      = reverse_dict(centered_dict,
                                           self.geometry_format_string[4])

    self.float_type         = reverse_dict(float_dict,
                                           self.geometry_format_string[5])

    self.has_iblank         = reverse_dict(iblank_dict,
                                           self.geometry_format_string[6])

    # *************************************************************************
    # read zone info (starting at index 1100)
    if self.is_structured == True:

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        if self.geometry_time_type == 'constant':

            for nz in range(self.n_zones):
                # instantiate zone and read info from file
                zone = StructuredZone()

                zone.geometry_header_length = \
                    structured_header_length[self.geometry_time_type]

                zone.number = len(self.zones)

                # reads geometry zone name (32 bytes)
                name = read_string(bytes_data,
                                   1100 + nz*zone.geometry_header_length, 32)
                zone._set_string(name, 'name', 32)
                zone._set_string(name, 'geometry_name', 32)

                # reads structured dimensions
                zone.iMax = read_int(bytes_data,
                                     1100 + 32 + nz*zone.geometry_header_length)
                zone.jMax = read_int(bytes_data,
                                     1100 + 36 + nz*zone.geometry_header_length)

                zone._update_geometry_info_str()
                self.zones.append(zone)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        elif self.geometry_time_type == 'aperiodic':

            for nz in range(self.n_zones):
                # instantiate zone and read info from file
                zone = StructuredZone()

                zone.geometry_header_length = \
                    structured_header_length[self.geometry_time_type]

                zone.number = len(self.zones)

                # reads geometry zone name (32 bytes)
                name = read_string(bytes_data,
                                   1100 + nz*zone.geometry_header_length, 32)
                zone._set_string(name, 'name', 32)
                zone._set_string(name, 'geometry_name', 32)

                # reads number of timesteps and structured dimensions
                zone.Nt = read_int(bytes_data,
                                   1100 + 32 + nz*zone.geometry_header_length)
                zone.iMax = read_int(bytes_data,
                                     1100 + 36 + nz*zone.geometry_header_length)
                zone.jMax = read_int(bytes_data,
                                     1100 + 40 + nz*zone.geometry_header_length)

                zone._update_geometry_info_str()
                self.zones.append(zone)

            # -----------------------------------------------------------------
            # Store 'Nt' in PWWPatch, check all zones have identical 'Nt'
            self.Nt = self.zones[0].Nt
            for z in range(self.n_zones):
                assert self.zones[z].Nt == self.Nt, \
                    "Error in file '{}': Zone {} has {} timesteps, while Zone 0 has {}!".format(geometry_filename, z, zone.Nt, self.Nt)
            # -----------------------------------------------------------------

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        elif self.geometry_time_type == 'periodic':

            for nz in range(self.n_zones):
                # instantiate zone and read info from file
                zone = StructuredZone()

                zone.geometry_header_length = \
                    structured_header_length[self.geometry_time_type]

                zone.number = len(self.zones)

                # reads geometry zone name (32 bytes)
                name = read_string(bytes_data,
                                   1100 + nz*zone.geometry_header_length, 32)
                zone._set_string(name, 'name', 32)
                zone._set_string(name, 'geometry_name', 32)

                # reads period, number of timesteps, and structured dimensions
                zone.period, _ = read_float(bytes_data,
                                            1100 + 32 + nz*zone.geometry_header_length)

                zone.Nt = read_int(bytes_data,
                                   1100 + 36 + nz*zone.geometry_header_length)

                zone.iMax = read_int(bytes_data,
                                     1100 + 40 + nz*zone.geometry_header_length)
                zone.jMax = read_int(bytes_data,
                                     1100 + 44 + nz*zone.geometry_header_length)

                zone._update_geometry_info_str()
                self.zones.append(zone)

            # -----------------------------------------------------------------
            # Store 'Nt' in PWWPatch, check all zones have identical 'Nt'
            self.Nt = self.zones[0].Nt
            for z in range(self.n_zones):
                assert self.zones[z].Nt == self.Nt, \
                    "Error in file '{}': Zone {} has {} timesteps, while Zone 0 has {}!".format(geometry_filename, z, zone.Nt, self.Nt)
            # -----------------------------------------------------------------
            # Store 'period' in PWWPatch, check all zones have identical 'period'
            self.period = self.zones[0].period
            for z in range(self.n_zones):
                assert self.zones[z].period == self.period, \
                    "Error in file '{}': Zone {} has period {} s, while Zone 0 has {}!".format(geometry_filename, z, zone.period, self.period)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        else:
            # TODO: implement other geometries header reader
            raise NotImplementedError("Can't read geometry header that is not constant, aperiodic or periodic - not implemented yet!")
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    # *************************************************************************
    else:
        # TODO: implement non-structured headers
        raise NotImplementedError("Can't read non-structured geometry zone info - not implemented yet!")

    # *************************************************************************
    assert (self.n_zones == len(self.zones)), \
        "Error in file '{}': Number of zones in format string doesn't match file data!".format(geometry_filename)

    # *************************************************************************


# #############################################################################
def _read_geometry_data(self, geometry_filename):
    """
    Reads PSU-WOPWOP geometry data from 'geometry_filename'.
    """

    # open fileand read binary content
    with open(geometry_filename, 'rb') as f:
        bytes_data = f.read()

    # *************************************************************************
    # If file is structured
    if self.is_structured == True:

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # constant geometry
        if self.geometry_time_type == 'constant':

            # start index for reading coordinates and normal vectors data
            field_start = 1100 + self.n_zones*self.zones[0].geometry_header_length

            # for each zone
            for nz in range(self.n_zones):

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # create empty numpy arrays for XYZ coordinates and normal
                # coordinates (and IBLANK data, if included)
                XYZ_coord = np.zeros((3, self.zones[nz].iMax, self.zones[nz].jMax),
                                     dtype=np.float32)

                normal_coord = np.zeros((3, self.zones[nz].iMax, self.zones[nz].jMax),
                             dtype=np.float32)

                self.zones[nz].add_StructuredConstantGeometry(XYZ_coord, normal_coord)

                if self.has_iblank == True:
                    self.zones[nz].geometry.iblank = \
                        np.zeros((self.zones[nz].iMax, self.zones[nz].jMax),
                                 dtype=np.int32)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # read XYZ coords and next index
                self.zones[nz].geometry.XYZ_coord, field_start = \
                    read_block(bytes_data, field_start, 3,
                               self.zones[nz].iMax, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # read normal vector coords and next index
                self.zones[nz].geometry.normal_coord, field_start = \
                    read_block(bytes_data, field_start, 3,
                               self.zones[nz].iMax, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # if file contains IBLANK (int) data, read that
                if self.has_iblank == True:
                    self.zones[nz].geometry.iblank, field_start = \
                        read_IBLANKblock(bytes_data, field_start,
                                         self.zones[nz].iMax,
                                         self.zones[nz].jMax)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # aperiodic geometry
        elif self.geometry_time_type == 'aperiodic':

            # start index for reading coordinates and normal vectors data
            field_start = 1100 + self.n_zones*self.zones[0].geometry_header_length

            # -----------------------------------------------------------------
            # for each timestep...
            for nt in range(self.Nt):

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # for each zone
                for nz in range(self.n_zones):

                    # .........................................................
                    # create empty numpy arrays for XYZ coordinates, normal
                    # coordinates, IBLANK data (if included), and time steps
                    XYZ_coord = np.zeros((self.Nt, 3, self.zones[nz].iMax, self.zones[nz].jMax),
                                         dtype=np.float32)

                    normal_coord = np.zeros((self.Nt, 3, self.zones[nz].iMax, self.zones[nz].jMax),
                                            dtype=np.float32)

                    self.zones[nz].add_StructuredAperiodicGeometry(XYZ_coord, normal_coord)

                    if self.has_iblank == True:
                        self.zones[nz].geometry.iblank = \
                            np.zeros((self.Nt, self.zones[nz].iMax, self.zones[nz].jMax),
                                     dtype=np.int32)

                    self.zones[nz].time_steps = np.zeros((self.Nt,), dtype=np.float32)

                    # .........................................................
                    # read current time value and next index
                    self.zones[nz].time_steps[nt], field_start = read_float(bytes_data, field_start)

                    # read XYZ coords and next index
                    self.zones[nz].geometry.XYZ_coord[nt, :, :, :], field_start = \
                        read_block(bytes_data, field_start, 3,
                                   self.zones[nz].iMax, self.zones[nz].jMax)

                    # read normal vector coords and next index
                    self.zones[nz].geometry.normal_coord[nt, :, :, :], field_start = \
                        read_block(bytes_data, field_start, 3,
                                   self.zones[nz].iMax, self.zones[nz].jMax)

                    # if file contains IBLANK (int) data, read that
                    if self.has_iblank == True:
                        self.zones[nz].geometry.iblank[nt, :, :, :], field_start = \
                            read_IBLANKblock(bytes_data, field_start,
                                             self.zones[nz].iMax,
                                             self.zones[nz].jMax)
                    # .........................................................

            # -----------------------------------------------------------------
            # compare time_steps across all zones, ensure they are identical
            self.time_steps = np.copy(self.zones[0].time_steps)
            for nz in range(self.n_zones):
                assert np.allclose(self.zones[nz].time_steps, self.time_steps),\
                    "Error reading file {}: Zone {} time steps do not match Zone 0 time steps!".format(geometry_filename, nz)
            # -----------------------------------------------------------------

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # periodic geometry
        elif self.geometry_time_type == 'periodic':

            # start index for reading coordinates and normal vectors data
            field_start = 1100 + self.n_zones*self.zones[0].geometry_header_length

            # -----------------------------------------------------------------
            # for each timestep...
            for nt in range(self.Nt):

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
                # for each zone
                for nz in range(self.n_zones):

                    # .........................................................
                    # create empty numpy arrays for XYZ coordinates, normal
                    # coordinates, IBLANK data (if included), and time steps
                    XYZ_coord = np.zeros((self.Nt, 3, self.zones[nz].iMax, self.zones[nz].jMax),
                                         dtype=np.float32)

                    normal_coord = np.zeros((self.Nt, 3, self.zones[nz].iMax, self.zones[nz].jMax),
                                            dtype=np.float32)

                    self.zones[nz].add_StructuredPeriodicGeometry(XYZ_coord, normal_coord,
                                                                  self.period)

                    if self.has_iblank == True:
                        self.zones[nz].geometry.iblank = \
                            np.zeros((self.Nt, self.zones[nz].iMax, self.zones[nz].jMax),
                                     dtype=np.int32)

                    self.zones[nz].time_steps = np.zeros((self.Nt,), dtype=np.float32)

                    # .........................................................
                    # read current time value and next index
                    self.zones[nz].time_steps[nt], field_start = read_float(bytes_data, field_start)

                    # read XYZ coords and next index
                    self.zones[nz].geometry.XYZ_coord[nt, :, :, :], field_start = \
                        read_block(bytes_data, field_start, 3,
                                   self.zones[nz].iMax, self.zones[nz].jMax)

                    # read normal vector coords and next index
                    self.zones[nz].geometry.normal_coord[nt, :, :, :], field_start = \
                        read_block(bytes_data, field_start, 3,
                                   self.zones[nz].iMax, self.zones[nz].jMax)

                    # if file contains IBLANK (int) data, read that
                    if self.has_iblank == True:
                        self.zones[nz].geometry.iblank[nt, :, :, :], field_start = \
                            read_IBLANKblock(bytes_data, field_start,
                                             self.zones[nz].iMax,
                                             self.zones[nz].jMax)
                    # .........................................................

            # -----------------------------------------------------------------
            # compare time_steps across all zones, ensure they are identical
            self.time_steps = np.copy(self.zones[0].time_steps)
            for nz in range(self.n_zones):
                assert np.allclose(self.zones[nz].time_steps, self.time_steps),\
                    "Error reading file {}: Zone {} time steps do not match Zone 0 time steps!".format(geometry_filename, nz)
            # -----------------------------------------------------------------

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        else:
            # TODO: implement other geometries data reader
            raise NotImplementedError("Can't read geometry data that is not constant, aperiodic or periodic - not implemented yet!")
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    # *************************************************************************
    else:
        # TODO: implement non-structured data reader
        raise NotImplementedError("Can't read non-structured geometry zone data - not implemented yet!")

    # *************************************************************************


# #############################################################################
# %%PSU-WOPWOP geometry file writers
# #############################################################################

def _write_geometry_header(self, geometry_filename):
    """
    Writes PSU-WOPWOP geometry header to 'geometry_filename'.
    """

    # build format string with most recent format values
    self._build_geometry_format_string()

    # open or create new file
    with open(geometry_filename, 'wb+') as file:

        # write 'magic number' 42 to first 4 bytes
        write_binary(file,MAGICNUMBER)

        # write version number
        write_binary(file,self.version_number_major)
        write_binary(file,self.version_number_minor)

        # write units string (32 chars)
        write_string(file,self.units_string, 32)

        # write comments string (1024 bytes)
        write_string(file,self.geometry_comment, 1024)

        # write format string (8 ints, 32 bytes)
        for n in range(8):
            write_binary(file,self.geometry_format_string[n])

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # write zone info
        if self.is_structured == True:

            # -----------------------------------------------------------------
            if self.geometry_time_type == 'constant':

                # for each zone...
                for nz in range(self.n_zones):

                    # write name (32-byte string)
                    write_string(file,self.zones[nz].geometry_name, 32)

                    # write iMax and jMax (4 byte ints)
                    write_binary(file,self.zones[nz].iMax)
                    write_binary(file,self.zones[nz].jMax)

            # -----------------------------------------------------------------
            elif self.geometry_time_type == 'aperiodic':

                # for each zone...
                for nz in range(self.n_zones):

                    # write name (32-byte string)
                    write_string(file,self.zones[nz].geometry_name, 32)

                    # write number of timesteps
                    write_binary(file,self.Nt)

                    # write iMax and jMax (4 byte ints)
                    write_binary(file,self.zones[nz].iMax)
                    write_binary(file,self.zones[nz].jMax)

            # -----------------------------------------------------------------
            elif self.geometry_time_type == 'periodic':

                # for each zone...
                for nz in range(self.n_zones):

                    # write name (32-byte string)
                    write_string(file, self.zones[nz].geometry_name, 32)

                    # write period
                    write_binary(file, self.period)

                    # write number of timesteps
                    write_binary(file, self.Nt)

                    # write iMax and jMax (4 byte ints)
                    write_binary(file, self.zones[nz].iMax)
                    write_binary(file, self.zones[nz].jMax)

            # -----------------------------------------------------------------
            else:
                # TODO: implement other geometries header writer
                raise NotImplementedError("Can't write geometry header that is not constant, aperiodic or periodic - not implemented yet!")
            # -----------------------------------------------------------------

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        else:
            # TODO: implement non-structured headers
            raise NotImplementedError("Can't write non-structured geometry zone info - not implemented yet!")
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


# #############################################################################
def _write_geometry_data(self, geometry_filename):
    """
    Writes PSU-WOPWOP geometry patch data to 'geometry_filename'
    Must be immediately preceded by a call to
    '_write_geometry_header(geometry_filename)'!
    """

    # open file in append mode - no need to adjust index
    with open(geometry_filename, 'ab') as file:

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # If file is structured
        if self.is_structured == True:

            # -----------------------------------------------------------------
            # constant geometry
            if self.geometry_time_type == 'constant':

                # for each zone
                for nz in range(self.n_zones):

                    # write geometry data (and iBlank data, if any)
                    write_block(file,self.zones[nz].geometry.XYZ_coord)
                    write_block(file,self.zones[nz].geometry.normal_coord)
                    if self.has_iblank == True:
                        write_block(file,self.zones[nz].geometry.iblank)

            # -----------------------------------------------------------------
            # aperiodic geometry
            elif self.geometry_time_type == 'aperiodic':

                # for each timestep...
                for nt in range(self.Nt):

                    # for each zone
                    for nz in range(self.n_zones):

                        # write current time step and geometry data
                        write_binary(file, self.time_steps[nt])
                        write_block(file, self.zones[nz].geometry.XYZ_coord)
                        write_block(file, self.zones[nz].geometry.normal_coord)

                        if self.has_iblank == True:
                            write_block(file,self.zones[nz].geometry.iblank)

            # -----------------------------------------------------------------
            # periodic geometry
            elif self.geometry_time_type == 'periodic':

                # for each timestep...
                for nt in range(self.Nt):

                    # for each zone
                    for nz in range(self.n_zones):

                        # write current time step and geometry data
                        write_binary(file, self.time_steps[nt])
                        write_block(file, self.zones[nz].geometry.XYZ_coord)
                        write_block(file, self.zones[nz].geometry.normal_coord)

                        if self.has_iblank == True:
                            write_block(file,self.zones[nz].geometry.iblank)

            # -----------------------------------------------------------------
            else:
                # TODO: implement other geometries data writer
                raise NotImplementedError("Can't write geometry data that is not constant, aperiodic or periodic - not implemented yet!")
            # -----------------------------------------------------------------

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        else:
            # TODO: implement non-structured headers
            raise NotImplementedError("Can't write non-structured geometry data - not implemented yet!")
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

