# -*- coding: utf-8 -*-
"""
Collection of convenience routines to parse and create PSU-WOPWOP input files
version 1.0.

Author:
    Fabio Casagrande Hirono
    Aug 2021
"""


import struct
import numpy as np


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# PySU-WOPWOP program constants and dicts

MAGICNUMBER = 42

# whether file is geometry file or node subset geometry file
geom_dict = {'geometry'     : 1,
             'node_subset'  :-1}

# whether file is structured or not
structured_dict = {'y' : 1,
                   'n' : 2}

# type of geometry (temporal behaviour)
# --> for functional files:
time_func_dict = {'constant'                 :1,
                  'periodic'                 :2,
                  'aperiodic'                :3,
                  'mult_time_aperiodic'      :4}

# for geometry files: all of functional types, plus quasiperiodic types
time_geom_dict = {**time_func_dict,
                  'quasiperiodic'            :5,
                  'mult_time_quasiperiodc'   :6}

# whether normal vectors and areas are node-centered or face-centered
centered_dict = {'node' : 1,
                 'face' : 2}

# what type of loading data is contained in functional file
data_dict = {'surf_pressure'    : 1,
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
iblank_dict = {'y' : 1,
               'n' : 0}

# Digit reserved for future use - must be '0' in this version
RESERVED_DIGIT = 0


# %% ############################################################################
class InputFile:
    """
    Parent class to read, write and store PSU-WOPWOP input files.

    This is subclassed to obtain separate classes for functional loading files
    and for geometry patch files.
    """

    def __init__(self):
        # PSU-WOPWOP input file v1.0 by default
        self.version_number1 = 1
        self.version_number2 = 0

        # values are little endian, 4-byte, signed by default
        self.endianness = 'little'
        self.value_length = 4
        self.is_signed = True

        # File type description (empty at initialization)
        self.n_zones = 0
        self.is_structured = ''
        self.time_type = ''
        self.centered_type = ''
        self.float_type = ''

        # 1024-byte comment string
        self.comment_string = ''

        # list of zones' information (stored as subclasses of Zone class)
        self.zones = []


    def _read_XYZblock(self, bytes_data, start_index, num_dims, iMax, jMax):
        """
        Reads a block of XYZ coordinates in PLOT3D format from a binary file.

        The block of data is a (num_dims, iMax, jMax)-shaped array of float32
        values.

        Parameters
        ----------
        bytes_data : file object
            File object obtained from calling 'open(filename)'.

        start_index : int
            Index indicating initial position of XYZ block within binary file.

        num_dims : int
            Number of dimensions to be read - e.g. '3' for XYZ, '2' for XY

        iMax : int
            First dimension of the data block

        jMax : int
            Second dimension of the data block

        Returns
        -------
        XYZ_data : (num_dims, iMax, jMax) array_like
            Numpy array containing the block data using float32 numbers.

        next_index : int
            Index for the next data field in the binary file.


        Notes
        -----
        The data in the file is assumed to be organized as:

        #####################################################################
        -Dim 0 (i.e X):

            X(i=0, j=0),    X(i=1, j=0),    ...     X(i=iMax, j=0)
            X(i=0, j=1),    ...             ...     X(i=iMax, j=1)
            ...                                     ...
            X(i=0, j=jMax)  ...                     X(i=iMax, j=jMax)
        #####################################################################
        -Dim 1 (i.e. Y):

            Y(i=0, j=0),    Y(i=1, j=0),    ...     Y(i=iMax, j=0)
            Y(i=1, j=0),    ...

            ...etc...

        #####################################################################
        """

        XYZ_data = np.zeros((num_dims, iMax, jMax), dtype=np.float32)

        # read XYZ coords
        for n in range(num_dims):
            for j in range(jMax):
                for i in range(iMax):

                    # fields are (i,j) order, hence j*iMax + i
                    current_index = (start_index
                                     + ((n*iMax*jMax + j*iMax + i)
                                        *self.value_length))

                    XYZ_data[n, i, j] = _read_float(bytes_data, current_index,
                                                    self.value_length, self.endianness)

        # increase start index by ('value_length' bytes * num_dims coords * iMax * jMax)
        next_index = start_index + self.value_length*(num_dims*iMax*jMax)

        return XYZ_data, next_index


# %% ############################################################################
class LoadingFunctionalFile(InputFile):
    """
    Child class to read, write and store PSU-WOPWOP loading file v1.0
    information.

    Note: Pressure input is the *gage* pressure, NOT the absolute pressure!
    """

    # **********************************************************************
    def __init__(self, filename=None):

        # initialize superclass
        super().__init__()

        self.data_type = ''
        self.ref_frame_type = ''

        self.zones_with_data = []
        self.n_zones_with_data = 0

        # if a file name is given, read that file to current instance
        if filename:
            self.filename = filename
            self.read_functional_file(filename)


    # **********************************************************************
    def _build_format_string(self):
        """
        Create list of data format values, later written to file as string of ints.
        """

        self.format_string = []
        self.format_string.append(2)               # indicate functional data file
        self.format_string.append(self.n_zones)
        self.format_string.append(structured_dict[self.is_structured])
        self.format_string.append(time_func_dict[self.time_type])
        self.format_string.append(centered_dict[self.centered_type])
        self.format_string.append(data_dict[self.data_type])
        self.format_string.append(ref_frame_dict[self.ref_frame_type])
        self.format_string.append(float_dict[self.float_type])
        self.format_string.append(RESERVED_DIGIT)   # reserved digit '0' inserted twice
        self.format_string.append(RESERVED_DIGIT)


    # **********************************************************************
    def print_info(self):
        """
        Prints a summary of the file info
        """

        if hasattr(self, 'filename'):
            print('File {} is a loading functional file for PSU-WOPWOP v{}.{}'.format(self.filename, self.version_number1, self.version_number2))
        else:
            print('File is a loading functional file for PSU-WOPWOP v{}.{}'.format(self.version_number1, self.version_number2))

        print('Comments: ', self.comment_string)

        print('Parsing format string:')
        print('\t--> Number of zones:               ', self.n_zones)
        print('\t--> Is file structured:            ', self.is_structured)
        print('\t--> Temporal behaviour:            ', self.time_type)
        print('\t--> Vectors/areas are centred at:  ', self.centered_type)
        print('\t--> Type of pressure/loading data: ', self.data_type)
        print('\t--> Reference frame:               ', self.ref_frame_type)
        print('\t--> Float data precision:          ', self.float_type)

        print('\n')
        print('Parsing zone specification:')

        print('No. zones with functional data: ', self.n_zones_with_data)
        for zone in self.zones:
            print('\t--> Zone name: ', zone.name)
            print('\t--> iMax: ', zone.iMax)
            print('\t--> jMax: ', zone.jMax)

        print('\n\n')


    # **********************************************************************
    def addStructuredConstantZone(self, name, data, zone_number):
        """
        Adds structured, constant zone containing pressure or loading vector data.

            name:           32-byte string

        If zone contains pressure info:
            data:       (iMax, jMax) array

        If zone contains loading vetors:
            data:       (3, iMax, jMax) array

            zone_number: int with zone index in geom file


        Parameters
        ----------
        name : string
            Zone name as a 32-byte string.

        data : (iMax, jMax) or (3, iMax, jMax) array_like
            The array of data to be added. Its shape is (iMax, jMax) for
            pressure data, and (3, iMax, jMax) for loading vector data.

        zone_number : int
            Zone number, as input to geometry patch file. For example, '5' if
            this zone is the 5th zone inserted to the patch file.

        Returns
        -------
        None.

        """

        # instantiate new zone
        zone = StructuredConstantZone(self.value_length)
        zone.name = name[:32]           # cap 'name' length to 32 bytes

        zone.number = self.n_zones+1

        # ------------------------------------------------------------------
        # if data is (iMax, jMax)-shaped  array of surface pressures
        if self.data_type == 'surf_pressure':
            zone.iMax = data.shape[0]
            zone.jMax = data.shape[1]

            # create empty numpy arrays for surface pressures, copy input data
            zone.pressures = np.zeros((zone.iMax, zone.jMax), dtype=np.float32)
            zone.pressures = np.copy(data)

            self.zones.append(zone)
            self.n_zones += 1
            self.n_zones_with_data += 1
            self.zones_with_data.append(zone_number)

        # ------------------------------------------------------------------
        # if data is (3, iMax, jMax)-shaped array of surface loading vectors
        elif self.data_type == 'surf_loading_vec':
            zone.iMax = data.shape[1]
            zone.jMax = data.shape[2]

            # create empty numpy arrays for surface loading vectors, copy input data
            zone.loading_vectors = np.zeros((3, zone.iMax, zone.jMax), dtype=np.float32)
            zone.loading_vectors = np.copy(data)

            self.zones.append(zone)
            self.n_zones += 1
            self.n_zones_with_data += 1

        # ------------------------------------------------------------------


    # **********************************************************************
    def write_functional_file(self, filename):
        """ Writes loading functional data to 'filename'."""
        self._write_header(filename)
        self._write_data(filename)

    def read_functional_file(self, filename):
        """ Reads 'filename' for loading functional data, and print a summary. """
        self._read_header(filename)
        self._read_data(filename)
        self.print_info()


    # **********************************************************************
    def _read_header(self, filename):
        """
        Reads PSU-WOPWOP loading functional header info from 'filename'.
        """

        # do initial check for 'magic number' and endianness
        self.endianness = initial_check(filename)

        # open fileand read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # check version number - must be '1'and '0'
        self.version_number1 = _read_int(bytes_data, 4, 4)
        self.version_number2 = _read_int(bytes_data, 8, 4)

        assert ((self.version_number1 == 1) and (self.version_number2 == 0)), 'File version is not v1.0!'

        # read comments string (1024 bytes, starting at index 12)
        self.comment_string = _read_string(bytes_data, 12, 1024)


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read format string (10 ints, 40 bytes, starting at index 1036)
        self.format_string = []
        for n in range(10):
            self.format_string.append(_read_int(bytes_data, 1036 + 4*n,
                                                self.value_length, self.endianness))

        # Populate file type description
        # --> format_string[0] is '2' to indicate functional data file
        assert (self.format_string[0] == 2), "Format string does not start with '2', not a PSU-WOPWOP functional file!"

        self.n_zones            = self.format_string[1]
        self.is_structured      = _reverse_dict(structured_dict, self.format_string[2])
        self.time_type          = _reverse_dict(time_func_dict, self.format_string[3])
        self.centered_type      = _reverse_dict(centered_dict, self.format_string[4])
        self.data_type          = _reverse_dict(data_dict, self.format_string[5])
        self.ref_frame_type     = _reverse_dict(ref_frame_dict, self.format_string[6])
        self.float_type         = _reverse_dict(float_dict, self.format_string[7])
        # --> format_string[8] is reserved for future use, and must be '0' in this version
        # --> format_string[9] is reserved for future use, and must be '0' in this version


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # zone specification

        # read number of zones containing functional data (4 bytes)
        self.n_zones_with_data = _read_int(bytes_data, 1076, 4)

        # create list of zones containing functional data
        # -->>  negative numbers indicate zones for which WOPWOP should *NOT*
        #       calculate thickness noise - e.g. loading patch
        self.zones_with_data = []
        for z in range(self.n_zones_with_data):
            self.zones_with_data.append(_read_int(bytes_data, 1080 + z*4, 4))


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read header info

        zone_info_start = 1080 + self.n_zones_with_data*4

        # version 1.0
        if self.version_number2 == 0:

            # ------------------------------------------------------------------
            # structured, constant geometry
            if self.is_structured == 'y':
                if self.time_type == 'constant':

                    for nz in range(self.n_zones_with_data):
                        # instantiate zone and read info from file
                        zone = StructuredConstantZone(self.value_length)

                        zone.name = _read_string(bytes_data, zone_info_start + nz*zone.header_length, 32)
                        zone.iMax = _read_int(bytes_data, zone_info_start + 32 + nz*zone.header_length, 4)
                        zone.jMax = _read_int(bytes_data, zone_info_start + 36 + nz*zone.header_length, 4)

                        self.zones.append(zone)

                else:
                    # ------------------------------------------------------------------
                    # TODO: implement non-constant loading
                    print('Reading non-constant functional zone info not implemented yet!')

            # ------------------------------------------------------------------
            else:
                # TODO: implement non-structured headers
                print('Reading non-structured functional zone info not implemented yet!')

        # ------------------------------------------------------------------
        # TODO: implement header reader for version 1.1
        elif self.version_number2 == '1':
            print('Reading functional data file v1.1 header not implemented yet!')

        # ------------------------------------------------------------------


    # **********************************************************************
    def _write_header(self, filename):
        """
        Writes PSU-WOPWOP functional file header to 'filename'
        """

        # build format string with most recent format values
        self._build_format_string()

        with open(filename, 'wb') as f:

            # write 'magic number' 42 to first 4 bytes
            _write_binary(f, MAGICNUMBER, self.value_length, self.endianness,
                          self.is_signed)

            # write version number
            _write_binary(f, self.version_number1, self.value_length, self.endianness,
                          self.is_signed)
            _write_binary(f, self.version_number2, self.value_length, self.endianness,
                          self.is_signed)

            # write comments string (1024 bytes)
            _write_string(f, self.comment_string, 1024)

            # write format string (10 ints, 40 bytes)
            for n in range(10):
                _write_binary(f, self.format_string[n], self.value_length,
                              self.endianness, self.is_signed)

            #****************************************************************
            # write zone specification

            # write number of zones with data
            _write_binary(f, self.n_zones_with_data, self.value_length,
                          self.endianness, self.is_signed)

            # write list of those zones
            for z in range(self.n_zones_with_data):
                _write_binary(f, self.zones_with_data[z], self.value_length,
                              self.endianness, self.is_signed)

            #****************************************************************
            # write zone info

            if self.is_structured == 'y':
                if self.time_type == 'constant':

                    # for each zone containing data...
                    for nz in range(self.n_zones_with_data):

                        # write name (32-byte string)
                        _write_string(f, self.zones[nz].name, 32)

                        # write iMax and jMax (4 byte ints)
                        _write_binary(f, self.zones[nz].iMax)
                        _write_binary(f, self.zones[nz].jMax)

                else:
                    # TODO: implement non-constant functional data header
                    print('Writing non-constant functional data header not implemented yet!')

            else:
                # TODO: implement non-structured headers
                print('Writing non-structured functional data header not implemented yet!')


    # **********************************************************************
    def _read_data(self, filename):
        """
        Reads PSU-WOPWOP loading functional data from 'filename'.
        """

        # open file and read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # structured, constant geometry
        if self.is_structured == 'y':
            if self.time_type == 'constant':

                # start index for reading functional data
                field_start = (1080 + self.n_zones_with_data*4
                               + self.n_zones_with_data*self.zones[0].header_length)

                # -----------------------------------------------------------
                # if data is surface pressure
                if self.data_type == 'surf_pressure':

                    # for each zone
                    for nz in range(self.n_zones_with_data):

                        # create empty numpy arrays for pressure data
                        self.zones[nz].pressures = np.zeros((self.zones[nz].iMax, self.zones[nz].jMax),
                                                            dtype=np.float32)

                        # read pressure data and next index
                        self.zones[nz].pressures, field_start = self._read_XYZblock(bytes_data, field_start,
                                                                                    1, self.zones[nz].iMax, self.zones[nz].jMax)

                # -----------------------------------------------------------
                # if data is loading vectors
                elif self.data_type == 'surf_loading_vec':

                    # for each zone
                    for nz in range(self.n_zones_with_data):

                        # create empty numpy arrays for pressure data
                        self.zones[nz].loading_vectors = np.zeros((3, self.zones[nz].iMax, self.zones[nz].jMax),
                                                            dtype=np.float32)

                        # read pressure data and next index
                        self.zones[nz].loading_vectors, field_start = self._read_XYZblock(bytes_data, field_start,
                                                                                          3, self.zones[nz].iMax, self.zones[nz].jMax)
                # -----------------------------------------------------------
                elif self.data_type == 'flow_params':
                    print('Reading flow parameters not implemented yet!')
                # -----------------------------------------------------------

            else:
                # TODO: read non-constant functional data
                print('Reading non-constant functional data not implemented yet!')
        else:
            # TODO: read non-structured functional data
            print('Reading non-structured functional data not implemented yet!')


    # **********************************************************************
    def _write_data(self, filename):
        """
        Writes PSU-WOPWOP functional data to 'filename'
        Must be immediately preceded by a call to '_write_header(filename)'!
        """

        # open file in append mode - no need to adjust index
        with open(filename, 'ab') as f:

            # -----------------------------------------------------------
            # constant geometry
            if self.time_type == 'constant':

                # -----------------------------------------------------------
                # if data is surface pressure
                if self.data_type == 'surf_pressure':

                    # for each zone
                    for nz in range(self.n_zones_with_data):

                        # write pressure data in Fortran (column-major) order
                        for j in range(self.zones[nz].jMax):
                            for i in range(self.zones[nz].iMax):
                                _write_binary(f, self.zones[nz].pressures[i, j])

                # -----------------------------------------------------------
                # if data is surface loading vectors
                elif self.data_type == 'surf_loading_vec':

                    # for each zone
                    for nz in range(self.n_zones_with_data):

                        # write loading vectors in Fortran (column-major) order
                        for n in range(3):
                            for j in range(self.zones[nz].jMax):
                                for i in range(self.zones[nz].iMax):
                                    _write_binary(f, self.zones[nz].loading_vectors[n, i, j])

                # -----------------------------------------------------------
                elif self.data_type == 'flow_params':
                    # TODO: write flow params data
                    print('Reading flow parameters not implemented yet!')

            # -----------------------------------------------------------
            else:
                # TODO: write non-constant functional data
                print('Writing non-constant functional data not implemented yet!')
            # -----------------------------------------------------------


# %% #########################################################################
class GeometryPatchFile(InputFile):
    """
    Child class to read, write and store PSU-WOPWOP geometry patch file v1.0
    information
    """

    def __init__(self, filename=None):

        # initialize superclass
        super().__init__()

        self.geom_type = ''
        self.iblank_included = ''

        # 32-byte string describing units of pressure to be used - e.g. 'Pa'
        self.units_string = ''

        # format string (8 ints)
        self.format_string = []

        # if a file name is given, read that file to current instance
        if filename:
            self.filename = filename
            self.read_patch_file(filename)

    # **********************************************************************
    def _build_format_string(self):
        """
        Create list of data format values, later written to file as string of ints.
        """

        self.format_string = []
        self.format_string.append(geom_dict[self.geom_type])
        self.format_string.append(self.n_zones)
        self.format_string.append(structured_dict[self.is_structured])
        self.format_string.append(time_geom_dict[self.time_type])
        self.format_string.append(centered_dict[self.centered_type])
        self.format_string.append(float_dict[self.float_type])
        self.format_string.append(iblank_dict[self.iblank_included])
        self.format_string.append(RESERVED_DIGIT)


    # **********************************************************************
    def print_info(self):
        """
        Prints a summary of the patch file info
        """

        if hasattr(self, 'filename'):
            print('File {} is a geometry patch file for PSU-WOPWOP v{}.{}'.format(self.filename, self.version_number1, self.version_number2))
        else:
            print('File is a geometry patch file for PSU-WOPWOP v{}.{}'.format(self.version_number1, self.version_number2))


        print('Units: ', self.units_string)

        print('Comments: ', self.comment_string)

        print('Parsing format string:')
        print('\t--> Geometry type:                 ', self.geom_type)
        print('\t--> Number of zones:               ', self.n_zones)
        print('\t--> Is file structured:            ', self.is_structured)
        print('\t--> Temporal behaviour:            ', self.time_type)
        print('\t--> Vectors/areas are centred at:  ', self.centered_type)
        print('\t--> Float data precision:          ', self.float_type)
        print('\t--> IBLANK values are included:    ', self.iblank_included)

        print('\n')
        print('Parsing zones info:')
        print('Number of zones: ', self.n_zones)

        for zone in self.zones:
            print('\n')
            print('\t--> Zone name: ', zone.name)
            print('\t--> iMax: ', zone.iMax)
            print('\t--> jMax: ', zone.jMax)

        print('\n\n')

    # **********************************************************************
    def addStructuredConstantZone(self, name, XYZ_coord, normal_coord):
        """
        Adds structured, constant zone containing mesh points and normal vectors.

        Parameters
        ----------
        name : string
            Zone name as a 32-byte string.

        XYZ_coord : (3, iMax, jMax) array_like
            Array of mesh point coordinates to be added.

        normal_coord : (3, iMax, jMax) array_like
            Array of normal vector coordinates to be added.

        zone_number : int
            Zone number, as input to geometry patch file. For example, '5' if
            this zone is the 5th zone inserted to the patch file.

        Returns
        -------
        None.

        """

        # instantiate new zone
        zone = StructuredConstantZone(self.value_length)

        # cap 'name' length to 32 bytes
        zone.name = name[:32]

        zone.number = self.n_zones+1
        zone.iMax = XYZ_coord.shape[1]
        zone.jMax = XYZ_coord.shape[2]

        # create empty numpy arrays for XYZ coordinates and normal
        # coordinates, copy input data
        zone.XYZ_coord    = np.zeros((3, zone.iMax, zone.jMax), dtype=np.float32)
        zone.normal_coord = np.zeros((3, zone.iMax, zone.jMax), dtype=np.float32)

        zone.XYZ_coord = np.copy(XYZ_coord)
        zone.normal_coord = np.copy(normal_coord)

        self.zones.append(zone)
        self.n_zones += 1


    # **********************************************************************
    def _read_IBLANKblock(self, bytes_data, start_index, iMax, jMax):
        """
        Reads a block of IBLANK data in PLOT3D format from a binary file.

        The block of data is a (iMax, jMax)-shaped array of int32 values.

        Parameters
        ----------
        bytes_data : file object
            File object obtained from calling 'open(filename)'.

        start_index : int
            Index indicating initial position of XYZ block within binary file.

        iMax : int
            First dimension of the data block

        jMax : int
            Second dimension of the data block

        Returns
        -------
        IBLANK_data : (iMax, jMax) array_like
            Numpy array containing IBLANK data as int32 numbers.

        next_index : int
            Index for the next data field in the binary file.

        Notes
        -----
        The data in the file is assumed to be organized as:

        #####################################################################
            X(i=0, j=0),    X(i=1, j=0),    ...     X(i=iMax, j=0)
            X(i=0, j=1),    ...             ...     X(i=iMax, j=1)
            ...                                     ...
            X(i=0, j=jMax)  ...                     X(i=iMax, j=jMax)
        #####################################################################

        """

        IBLANK_data = np.zeros((iMax, jMax), dtype=np.int32)

        for j in range(jMax):
            for i in range(iMax):

                # fields are (i,j) order, hence j*iMax + i
                current_index = (start_index
                                  + (j*iMax + i)*self.value_length)

                IBLANK_data[i, j] = _read_int(bytes_data, current_index,
                                              self.value_length, self.endianness)

        # increase start index by ('value_length' bytes * iMax * jMax)
        next_index = start_index + self.value_length*(iMax*jMax)

        return IBLANK_data, next_index


    # **********************************************************************
    def write_patch_file(self, filename):
        """ Writes geometry patch data to 'filename'."""
        self._write_header(filename)
        self._write_data(filename)


    def read_patch_file(self, filename):
        """ Reads 'filename' for geometry patch data, and prints a summary. """
        self._read_header(filename)
        self._read_data(filename)
        self.print_info()


    # **********************************************************************
    def _read_header(self, filename):
        """
        Reads PSU-WOPWOP geometry patch header info from 'filename'.
        """

        # do initial check for 'magic number' and endianness
        self.endianness = initial_check(filename)

        # open fileand read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # check version number - must be '1'and '0'
        self.version_number1 = _read_int(bytes_data, 4, 4)
        self.version_number2 = _read_int(bytes_data, 8, 4)

        assert ((self.version_number1 == 1) and (self.version_number2 == 0)), 'File version is not v1.0!'

        # read units string (32 chars, starting at index 12)
        self.units_string = _read_string(bytes_data, 12, 32)

        # read comments string (1024 bytes, starting at index 44)
        self.comment_string = _read_string(bytes_data, 44, 1024)

        # read format string (8 ints, 32 bytes, starting at index 1068)
        self.format_string = []
        for n in range(8):
            self.format_string.append(_read_int(bytes_data, 1068 + 4*n,
                                                self.value_length, self.endianness))

        # Populate file type description
        self.geom_type          = _reverse_dict(geom_dict, self.format_string[0])
        self.n_zones            = self.format_string[1]
        self.is_structured      = _reverse_dict(structured_dict, self.format_string[2])
        self.time_type          = _reverse_dict(time_geom_dict, self.format_string[3])
        self.centered_type      = _reverse_dict(centered_dict, self.format_string[4])
        self.float_type         = _reverse_dict(float_dict, self.format_string[5])
        self.iblank_included    = _reverse_dict(iblank_dict, self.format_string[6])


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read zone info (starting at index 1100)

        if self.is_structured == 'y':

            if self.time_type == 'constant':

                for nz in range(self.n_zones):
                    # instantiate zone and read info from file
                    zone = StructuredConstantZone(self.value_length)

                    zone.name = _read_string(bytes_data, 1100 + nz*zone.header_length, 32)
                    zone.iMax = _read_int(bytes_data, 1100 + 32 + nz*zone.header_length, 4)
                    zone.jMax = _read_int(bytes_data, 1100 + 36 + nz*zone.header_length, 4)

                    self.zones.append(zone)

            else:
                # TODO: implement non-constant geometries
                print('Reading non-constant geometry zone info not implemented yet!')

        else:
            # TODO: implement non-structured headers
            print('Reading non-structured geometry zone info not implemented yet!')


    # **********************************************************************
    def _write_header(self, filename):
        """
        Writes PSU-WOPWOP geometry patch header to 'filename'.
        """

        # build format string with most recent format values
        self._build_format_string()

        with open(filename, 'wb') as f:

            # write 'magic number' 42 to first 4 bytes
            _write_binary(f, MAGICNUMBER, self.value_length, self.endianness,
                          self.is_signed)

            # write version number
            _write_binary(f, self.version_number1, self.value_length, self.endianness,
                          self.is_signed)
            _write_binary(f, self.version_number2, self.value_length, self.endianness,
                          self.is_signed)

            # write units string (32 chars)
            _write_string(f, self.units_string, 32)

            # write comments string (1024 bytes)
            _write_string(f, self.comment_string, 1024)

            # write format string (8 ints, 32 bytes)
            for n in range(8):
                _write_binary(f, self.format_string[n], self.value_length,
                              self.endianness, self.is_signed)

            #****************************************************************
            # write zone info

            if self.is_structured == 'y':

                if self.time_type == 'constant':

                    # for each zone...
                    for nz in range(self.n_zones):

                        # write name (32-byte string)
                        _write_string(f, self.zones[nz].name, 32)

                        # write iMax and jMax (4 byte ints)
                        _write_binary(f, self.zones[nz].iMax)
                        _write_binary(f, self.zones[nz].jMax)

                else:
                    # TODO: implement non-constant geometries
                    print('Writing non-constant geometry zone info not implemented yet!')

            else:
                # TODO: implement non-structured headers
                print('Writing non-structured geometry zone info not implemented yet!')


    # **********************************************************************
    def _read_data(self, filename):
        """
        Reads PSU-WOPWOP geometry patch data from 'filename'.
        """

        # open fileand read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # constant geometry
        if self.time_type == 'constant':

            # start index for reading coordinates and normal vectors data
            field_start = 1100 + self.n_zones*self.zones[0].header_length

            # ***************************************************************
            # for each zone
            for nz in range(self.n_zones):

                # -----------------------------------------------------------
                # create empty numpy arrays for XYZ coordinates and normal
                # coordinates (and IBLANK data, if included)
                self.zones[nz].XYZ_coord    = np.zeros((3, self.zones[nz].iMax, self.zones[nz].jMax), dtype=np.float32)
                self.zones[nz].normal_coord = np.zeros((3, self.zones[nz].iMax, self.zones[nz].jMax), dtype=np.float32)

                if self.iblank_included == 'y':
                    self.zones[nz].iblank = np.zeros((self.zones[nz].iMax, self.zones[nz].jMax), dtype=np.int32)

                # -----------------------------------------------------------
                # read XYZ coords and next index
                self.zones[nz].XYZ_coord, field_start = self._read_XYZblock(bytes_data, field_start,
                                                                            3, self.zones[nz].iMax, self.zones[nz].jMax)
                # -----------------------------------------------------------
                # read normal vector coords and next index
                self.zones[nz].normal_coord, field_start = self._read_XYZblock(bytes_data, field_start,
                                                                               3, self.zones[nz].iMax, self.zones[nz].jMax)
                # -----------------------------------------------------------
                # if file contains IBLANK (int) data, read that
                if self.iblank_included == 'y':
                    self.zones[nz].iblank, field_start = self._read_IBLANKblock(bytes_data, field_start,
                                                                                self.zones[nz].iMax, self.zones[nz].jMax)
                # -----------------------------------------------------------

        else:
            # TODO: read non-constant geometry data
            print('Reading non-constant geometry data not implemented yet!')


    # **********************************************************************
    def _write_data(self, filename):
        """
        Writes PSU-WOPWOP geometry patch data to 'filename'
        Must be immediately preceded by a call to '_write_header(filename)'!
        """

        # open file in append mode - no need to adjust index
        with open(filename, 'ab') as f:

            # -----------------------------------------------------------
            # constant geometry
            if self.time_type == 'constant':

                # for each zone
                for nz in range(self.n_zones):

                    #write XYZ coords
                    for n in range(3):
                        # write order is Fortran (column-major)
                        for j in range(self.zones[nz].jMax):
                            for i in range(self.zones[nz].iMax):
                                _write_binary(f, self.zones[nz].XYZ_coord[n, i, j])

                    # write normal vector coords
                    for n in range(3):
                        for j in range(self.zones[nz].jMax):
                            for i in range(self.zones[nz].iMax):
                                _write_binary(f, self.zones[nz].normal_coord[n, i, j])

                    # write IBLANK data
                    if self.iblank_included == 'y':
                        for j in range(self.zones[nz].jMax):
                            for i in range(self.zones[nz].iMax):
                                _write_binary(f, self.zones[nz].iblank[i, j])

            # -----------------------------------------------------------
            else:
                # TODO: write non-constant geometry data
                print('Writing non-constant geometry data not implemented yet!')
            # -----------------------------------------------------------


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

class Zone:
    """
    Parent class for zone data containing zone name, value length (typically 4
    bytes) and header length (fixed at 32 bytes)
    """
    def __init__(self, value_length):
        self.name = ''
        self.number = 0                     # zone numbers use 1-based index
        self.value_length = value_length    # length in bytes of values
        self.header_length = 32             # 'name' string has 32 bytes


# **********************************************************************
class StructuredZone(Zone):
    """
    Parent class for structured zone, containing structured dimensions (2 ints)
    """

    def __init__(self, value_length):
        super().__init__(value_length)
        self.iMax = 0
        self.jMax = 0
        self.header_length += 2*self.value_length


class PeriodicZone(Zone):
    """
    Parent class for time-varying period zone, containing period (in seconds)
    and number of time steps (integer)
    """
    def __init__(self, value_length):
        super().__init__(value_length)
        self.period = 0.
        self.nt = 0
        self.header_length += 2*self.value_length


class AperiodicZone(Zone):
    """
    Parent class for time-varying aperiod zone, containing number of time
    steps (integer)
    """
    def __init__(self, value_length):
        super().__init__(value_length)
        self.nt = 0
        self.header_length += value_length


# **********************************************************************
class StructuredConstantZone(StructuredZone):
    """
    Class for structured, constant zone - contains no time information
    """
    def __init__(self, value_length):
        super().__init__(value_length)


# TODO: add period and nt info on initialization
class StructuredPeriodicZone(StructuredZone, PeriodicZone):
    """
    Class for structured, periodic zone
    """
    def __init__(self, value_length):
        super().__init__(value_length)


# TODO: add nt info on initialization
class StructuredAperiodicZone(StructuredZone, AperiodicZone):
    """
    Class for structured, aperiodic zone
    """
    def __init__(self, value_length):
        super().__init__(value_length)



# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
# Assorted functions and checks

def _reverse_dict(my_dict, my_value):
    """
    Tries to find which key corresponds to a given value in a dictionary - i.e.
    reverses the dict lookup.

    This code assumes no two keys share the same value - as is the case with
    PSU-WOPWOP dicts used here.
    """

    return {value: key for key, value in my_dict.items()}[my_value]


def initial_check(filename):
    """
    Check the first 4 bytes of a file for the 'magic number' and return the
    file endianness. If the 'magic number' is not found, the file is probably
    not a PSU-WOPWOP file and an error is raised.
    """

    endianness_flag = 'little'

    # read first four bytes to check for 'magic number' 42 and endianness
    with open(filename, 'rb') as file:
        bytes_data = file.read(4)

    # if little endian, continue
    if bytes_data == b'*\x00\x00\x00':
        print('Magic number is correct - file is little endian\n')

    # if big endian, change flag and continue
    elif bytes_data == b'\x00\x00\x00*':
        endianness_flag = 'big'
        print('Magic number is correct - file is big endian\n')

    # if none, magic number is incorrect - probably not PSU-WOPWOP file!
    else:
        raise ValueError('Magic number is incorrect - file {} is probably not a PSU-WOPWOP patch file v1.0!'.format(filename))

    return endianness_flag


def _read_int(obj_name, start_index, n_bytes, endianness_flag='little'):
    """
    Reads one integer value from an open file and returns the unpacked value.

    Parameters
    ----------
    obj_name : bytes_object
        Object containing binary data, such as an open file object.

    start_index : int
        Starting index of value to be read.

    n_bytes : {1, 2, 4, 8}
        Size of the integer to be read, in bytes.

    endianness_flag : {'little', 'big'}, optional
        String indicating the byte endinaness to be used. The default is
        'little'.

    Returns
    -------
    out : int
        Integer value unpacked from file data.

    """

    n_bytes_dict = {1:'b', 2:'h', 4:'i', 8:'q'}

    endianness_dict = {'little':'<', 'big':'>'}

    return struct.unpack(endianness_dict[endianness_flag] + n_bytes_dict[n_bytes],
                         obj_name[start_index:start_index + n_bytes])[0]


def _read_float(obj_name, start_index, n_bytes, endianness_flag='little'):
    """
    Reads one float value from an open file and returns the unpacked value.

    Parameters
    ----------
    obj_name : bytes_object
        Object containing binary data, such as an open file object.

    start_index : int
        Starting index of value to be read.

    n_bytes : {4, 8}
        Size of the float to be read, in bytes.

    endianness_flag : {'little', 'big'}, optional
        String indicating the byte endinaness to be used. The default is
        'little'.

    Returns
    -------
    out : float
        Float value unpacked from file data.

    """

    n_bytes_dict = {4:'f', 8:'d'}

    endianness_dict = {'little':'<', 'big':'>'}

    return struct.unpack(endianness_dict[endianness_flag] + n_bytes_dict[n_bytes],
                         obj_name[start_index:start_index + n_bytes])[0]


def _write_binary(file, data, length=4, endianness_flag='little', is_signed=True):
    """
    Writes one value of data to an open binary file.

    Parameters
    ----------
    file : file object
        File object obtained from calling 'open(filename)'.

    data : int or float
        Value to be written to file. Must be int or float.

    length : {1, 2, 4, 8}, optional
        Byte length of the value to be written. The default is 4 bytes.

    endianness_flag : {'little', 'big'}, optional
        String indicating the byte endinaness to be used. The default is
        'little'.

    is_signed : boolean, optional
        Flag indicating whether the values are signed (True) or unsigned
        (False). The default is True.

    Returns
    -------
    None.

    Notes
    -----

    Data value to be written must be integer or floating point.

    Floating point data can only accept lengths of 4 and 8 bytes.

    """

    if type(data) is int:
        file.write(data.to_bytes(length, endianness_flag, signed=is_signed))

    # if data is python float or numpy float, write to file as float
    elif (isinstance(data, (float, np.floating))):

        endianness = {'little':'<', 'big':'>'}
        floatlen = {4:'f', 8:'d'}
        format_string = endianness[endianness_flag] + floatlen[length]

        file.write(struct.pack(format_string, data))


def _write_string(file, string, max_length):
    """
    Writes a ASCII-compatible string to an open binary file object, up to a
    maximum length. If string is shorter than 'max_length', pad with spaces.
    """

    # check string is ASCII compatible
    ascii_error = 'String is not ASCII compatible!'
    assert string[:max_length].isascii(), ascii_error

    # check string has length 'max_length', pad with spaces otherwise
    if len(string) < max_length:
        string += (max_length-len(string))*' '

    file.write(string[:max_length].encode('ascii'))


def _read_string(obj_name, start_index, len_string):
    """
    Reads strings of arbitrary length stored in binary format.
    """
    mystring = ''
    for i in range(len_string):
        mystring += chr(_read_int(obj_name, start_index + i, 1))

    return mystring

