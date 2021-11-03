# -*- coding: utf-8 -*-
"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

Author:
    Fabio Casagrande Hirono
    Oct 2021
"""


import struct
import numpy as np

# ##########################################################################
    # **********************************************************************
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
            # --------------------------------------------------------------
                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                    # ......................................................



# %% ##########################################################################
# PSU-WOPWOP program constants and dicts

MAGICNUMBER = 42

# values are little endian, 4-byte, signed by default
ENDIANNESS = 'little'
VALUE_LENGTH = 4
IS_SIGNED = True

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

# Digit reserved for future use - must be '0' in this version
RESERVED_DIGIT = 0


# %% ##########################################################################
class PWWPatch:
    """
    Parent class to read, write and store PSU-WOPWOP data and  files.
    """

    def __init__(self):
        # PSU-WOPWOP input file v1.0 by default
        self.version_number_major = 1
        self.version_number_minor = 0

        # File type description
        self.is_structured = True
        self.has_iblank = False
        self.float_type = 'single'
        self.centered_type = ''

        # 32-byte string describing units of pressure to be used - e.g. 'Pa'
        self.units_string = 'Pa'
        self.geometry_type = 'geometry'

        self.geometry_time_type = ''
        self.loading_time_type = ''

        self.loading_data_type = ''
        self.loading_ref_frame = ''

        # 1024-byte comment string for geometry file
        self.geometry_comment = ''

        # 1024-byte comment string for loading file
        self.loading_comment = ''

        self.geometry_format_string = []
        self.loading_format_string = []

        # list of all zones
        self.zones = []

        # list of zones containing loading data
        # --> use indices from self.zones
        # --> negative values indicate 'calc_thickness_noise' is False
        self.zones_with_loading_data = []

        self.n_zones = 0
        self.n_zones_with_loading_data = 0


    # **********************************************************************
    def print_info(self):
        """
        Prints a summary of the file info
        """

        print('\nPSU-WOPWOP version {}.{}'.format(self.version_number_major,
                                                  self.version_number_minor))

        print('\nGeneral info:')
        print('\t--> Number of zones:               ', len(self.zones))
        print('\t--> Is file structured:            ', self.is_structured)
        print('\t--> Vectors/areas are centred at:  ', self.centered_type)
        print('\t--> Float data precision:          ', self.float_type)

        print('\nGeometry info:')
        print('\t--> Geometry type:                 ', self.geometry_type)
        print('\t--> Geometry temporal behaviour:   ', self.geometry_time_type)
        print('\t--> IBLANK values are included:    ', self.has_iblank)
        print('\t--> Geometry comment string:\n\t"' +   self.geometry_comment + '"\n')

        print('\nLoading info:')
        print('\t--> Type of loading data:          ', self.loading_data_type)
        print('\t--> Loading temporal behaviour:    ', self.loading_time_type)
        print('\t--> Loading data reference frame:  ', self.loading_ref_frame)
        print('\t--> Loading comment string:\n\t"' +   self.loading_comment + '"\n')

        print('\nZone specification:')
        print('\tNo. zones:                 ', len(self.zones))
        print('\tNo. zones w/ loading data: ', len(self.zones_with_loading_data))

        for zone in self.zones:
            print(zone)
        print('\n\n')


    # **********************************************************************
    def write_geometry_file(self, geometry_filename):
        self._write_geometry_header(geometry_filename)
        self._write_geometry_data(geometry_filename)


    def write_loading_file(self, loading_filename):
        self._write_loading_header(loading_filename)
        self._write_loading_data(loading_filename)


    def read_geometry_file(self, geometry_filename):
        self._read_geometry_header(geometry_filename)
        self._read_geometry_data(geometry_filename)


    def read_loading_file(self, loading_filename):
        self._read_loading_header(loading_filename)
        self._read_loading_data(loading_filename)


    # **********************************************************************
    def add_StructuredZone(self, name, XYZ_coord, normal_coord,
                           calc_thickness_noise=True, loading_data=None):
        """
        Adds a new structured zone to PWWPatch instance. Expected argument
        array shapes and types will depend on attributes of current PWWPatch
        instance; see Notes below.

        Parameters
        ----------
        name : string
            32-byte long string containing zone name. It is used in both
            geometry and loading data fields by default.

        XYZ_coord : array_like

        normal_coord : array_like
            Array of mesh point coordinates. Its expected shape will
            depend on 'geometry_time_type' attribute of the parent PWWPatch
            instance; see Notes below.

        normal_coord : array_like
            Array of normal vector coordinates. Its expected shape will
            depend on 'geometry_time_type' attribute of the parent PWWPatch
            instance; see Notes below.

        calc_thickness_noise : boolean, optional
            Boolean flag denoting whether PSU-WOPWOP should calculate thickness
            noise generated by this zone. Default is True.

        loading_data : tuple (loading_data), optional
            Array containing loading data. Its expected shape will depend on
            'loading_time_type' and 'loading_data_type' attributes of the
            parent PWWPatch instance; see Notes below. Default is None.

        Returns
        -------
        None.

        Notes
        -----
        The type of input data for the geometry and loading will depend on the
        'geometry_time_type' and 'loading_time_type' attributes of the parent
        PWWPatch instance containing the zones. For the moment, only 'constant'
        data types are implemented.

        If  'geometry_time_type' is 'constant', the geometry arrays are:

            XYZ_coord : (3, iMax, jMax) array_like
                Array of mesh point coordinates to be added.

            normal_coord : (3, iMax, jMax) array_like
                Array of normal vector coordinates to be added.


        If 'loading time_type' is 'constant', the loading information is:

            loading_data : (iMax, jMax) or (3, iMax, jMax) array_like
                Array of constant loading data to be added. Its shape is (iMax,
                jMax) for pressure data, and (3, iMax, jMax) for loading vector
                data.

        """

        # check PWWPatch instance is indeed structured
        assert (self.is_structured == True), \
            "Cannot add structured zone - PWWPatch instance is not structured!"

        # instantiate new zone
        zone = StructuredZone()

        zone.iMax, zone.jMax = XYZ_coord.shape[1:]

        # cap 'name' length to 32 bytes
        zone.geometry_name = name[:32]
        zone.loading_name = name[:32]

        # zone number will always be the current length of 'zones'
        zone.number = len(self.zones)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # adds geometry data to zone
        if self.geometry_time_type == 'constant':
            zone.add_StructuredConstantGeometry(XYZ_coord, normal_coord)

        elif self.geometry_time_type == 'periodic':
            # TODO: implement structured periodic geometry
            print('Structured periodic geometry not implemented yet!')
            #zone.add_StructuredPeriodicGeometry(XYZ_coord, normal_coord)

        elif self.geometry_time_type == 'aperiodic':
            # TODO: implement structured aperiodic geometry
            print('Structured aperiodic geometry not implemented yet!')
            # zone.add_StructuredAperiodicGeometry(XYZ_coord, normal_coord)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # adds loading data, if there is any
        if loading_data is not None:

            # --------------------------------------------------------------
            if self.loading_time_type == 'constant':

                # check instance loading data type string vs. loading data array shape
                if self.loading_data_type == 'surf_pressure':
                    assert loading_data.shape == (zone.iMax, zone.jMax), \
                        "'loading_data' does not match expected shape for 'surf_pressure' (iMax, jMax)!"

                elif self.loading_data_type == 'surf_loading_vec':
                    assert loading_data.shape == (3, zone.iMax, zone.jMax), \
                        "'loading_data' does not match expected shape for 'surf_loading_vec' (3, iMax, jMax)!"

                elif self.loading_data_type == 'flow_params':
                    print("Can't check loading data shape for flow parameter data - not implemented yet!")

                zone.add_StructuredConstantLoading(loading_data, self.loading_data_type)

            # ----------------------------------------------------------------
            elif self.loading_time_type == 'periodic':

                # TODO: implement structured periodic loading
                print('Structured periodic loading not implemented yet!')
                #zone.add_StructuredPeriodicLoading(loading_data, self.loading_data_type)

            # ----------------------------------------------------------------
            elif self.loading_time_type == 'aperiodic':

                # TODO: implement structured aperiodic loading
                print('Structured aperiodic loading not implemented yet!')
                #zone.add_StructuredAperiodicLoading(loading_data, self.loading_data_type)

            # ----------------------------------------------------------------
            # set zone loading data flag
            zone.has_loading_data = True

            # set thickness noise flag, append zone number to
            # 'zones_with_loading_data'
            if calc_thickness_noise:
                zone.calc_thickness_noise = True
                self.zones_with_loading_data.append(zone.number)
            else:
                zone.calc_thickness_noise = False
                self.zones_with_loading_data.append(-zone.number)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # append new zone to zones list
        self.zones.append(zone)

        # update number of zones in PWWPatch instance
        self._update_n_zones()


    # **********************************************************************
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
            Number of dimensions to be read - e.g. '3' for XYZ, '2' for XY, or
            '1' for surface pressure data

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

        # read surface pressures - i.e. (iMax, jMax)-shaped
        if num_dims == 1:
            XYZ_data = np.zeros((iMax, jMax), dtype=np.float32)

            for j in range(jMax):
                for i in range(iMax):

                    # fields are (i,j) order, hence j*iMax + i
                    current_index = (start_index
                                     + ((j*iMax + i)*VALUE_LENGTH))

                    XYZ_data[i, j] = _read_float(bytes_data, current_index)

        # read XYZ or XY - i.e. (num_dims, iMax, jMax)-shaped
        else:
            XYZ_data = np.zeros((num_dims, iMax, jMax), dtype=np.float32)

            for n in range(num_dims):
                for j in range(jMax):
                    for i in range(iMax):

                        # fields are (i,j) order, hence j*iMax + i
                        current_index = (start_index
                                         + ((n*iMax*jMax + j*iMax + i)
                                            * VALUE_LENGTH))

                        XYZ_data[n, i, j] = _read_float(bytes_data, current_index)

        # increase start index by ('value_length' bytes * num_dims coords * iMax * jMax)
        next_index = start_index + VALUE_LENGTH*(num_dims*iMax*jMax)

        return XYZ_data, next_index


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

                IBLANK_data[i, j] = _read_int(bytes_data, current_index)

        # increase start index by ('value_length' bytes * iMax * jMax)
        next_index = start_index + VALUE_LENGTH*(iMax*jMax)

        return IBLANK_data, next_index


    # **********************************************************************
    def _read_geometry_header(self, filename):
        """
        Reads PSU-WOPWOP geometry header info from 'filename'.
        """

        # do initial check for 'magic number' and endianness
        file_endianness = initial_check(filename)
        assert file_endianness == ENDIANNESS, \
            "File endianness does not match pywopwop default!"

        # open fileand read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # check version number - must be '1'and '0'
        self.version_number_major = _read_int(bytes_data, 4)
        self.version_number_minor = _read_int(bytes_data, 8)

        assert ((self.version_number_major == 1)
                and (self.version_number_minor == 0)), \
            'File version is not v1.0!'

        # read units string (32 chars, starting at index 12)
        self.units_string = _read_string(bytes_data, 12, 32)

        # read comments string (1024 bytes, starting at index 44)
        self.geometry_comment = _read_string(bytes_data, 44, 1024)

        # read format string (8 ints, 32 bytes, starting at index 1068)
        self.geometry_format_string = []
        for n in range(8):
            self.geometry_format_string.append(_read_int(bytes_data, 1068 + 4*n))

        # Populate file type description
        self.geometry_type      = _reverse_dict(geom_dict, self.geometry_format_string[0])
        self.n_zones            = self.geometry_format_string[1]
        self.is_structured      = _reverse_dict(structured_dict, self.geometry_format_string[2])
        self.geometry_time_type = _reverse_dict(geometry_time_dict, self.geometry_format_string[3])
        self.centered_type      = _reverse_dict(centered_dict, self.geometry_format_string[4])
        self.float_type         = _reverse_dict(float_dict, self.geometry_format_string[5])
        self.has_iblank         = _reverse_dict(iblank_dict, self.geometry_format_string[6])


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read zone info (starting at index 1100)

        if self.is_structured == True:

            if self.geometry_time_type == 'constant':

                for nz in range(self.n_zones):
                    # instantiate zone and read info from file
                    zone = StructuredZone()

                    zone.geometry_name = _read_string(bytes_data, 1100 + nz*zone.header_length, 32)
                    zone.loading_name = ''
                    zone.iMax = _read_int(bytes_data, 1100 + 32 + nz*zone.header_length)
                    zone.jMax = _read_int(bytes_data, 1100 + 36 + nz*zone.header_length)

                    zone._update_geometry_info()
                    self.zones.append(zone)

            else:
                # TODO: implement non-constant geometries
                print('Reading non-constant geometry zone info not implemented yet!')

        else:
            # TODO: implement non-structured headers
            print('Reading non-structured geometry zone info not implemented yet!')

        assert (self.n_zones == len(self.zones)), \
            "Number of zones in format string doesn't match file data!"


    # **********************************************************************
    def _write_geometry_header(self, filename):
        """
        Writes PSU-WOPWOP geometry header to 'filename'.
        """

        # build format string with most recent format values
        self._build_geometry_format_string()

        with open(filename, 'wb') as f:

            # write 'magic number' 42 to first 4 bytes
            _write_binary(f, MAGICNUMBER)

            # write version number
            _write_binary(f, self.version_number_major)
            _write_binary(f, self.version_number_minor)

            # write units string (32 chars)
            _write_string(f, self.units_string, 32)

            # write comments string (1024 bytes)
            _write_string(f, self.geometry_comment, 1024)

            # write format string (8 ints, 32 bytes)
            for n in range(8):
                _write_binary(f, self.geometry_format_string[n])

            # --------------------------------------------------------------
            # write zone info
            if self.is_structured == True:

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                if self.geometry_time_type == 'constant':

                    # for each zone...
                    for nz in range(self.n_zones):

                        # write name (32-byte string)
                        _write_string(f, self.zones[nz].geometry_name, 32)

                        # write iMax and jMax (4 byte ints)
                        _write_binary(f, self.zones[nz].iMax)
                        _write_binary(f, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                else:
                    # TODO: implement non-constant geometries
                    print('Writing non-constant geometry zone info not implemented yet!')

            # --------------------------------------------------------------
            else:
                # TODO: implement non-structured headers
                print('Writing non-structured geometry zone info not implemented yet!')


    # **********************************************************************
    def _read_geometry_data(self, filename):
        """
        Reads PSU-WOPWOP geometry data from 'filename'.
        """

        # open fileand read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # constant geometry
        if self.geometry_time_type == 'constant':

            # start index for reading coordinates and normal vectors data
            field_start = 1100 + self.n_zones*self.zones[0].header_length

            # for each zone
            for nz in range(self.n_zones):

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                # create empty numpy arrays for XYZ coordinates and normal
                # coordinates (and IBLANK data, if included)
                XYZ_coord = np.zeros((3, self.zones[nz].iMax, self.zones[nz].jMax),
                                     dtype=np.float32)

                normal_coord = np.zeros((3, self.zones[nz].iMax, self.zones[nz].jMax),
                             dtype=np.float32)

                self.zones[nz].geometry = StructuredConstantGeometry(XYZ_coord, normal_coord)

                if self.has_iblank == True:
                    self.zones[nz].geometry.iblank = \
                        np.zeros((self.zones[nz].iMax, self.zones[nz].jMax),
                                 dtype=np.int32)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                # read XYZ coords and next index
                self.zones[nz].geometry.XYZ_coord, field_start = \
                    self._read_XYZblock(bytes_data, field_start, 3,
                                        self.zones[nz].iMax, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                # read normal vector coords and next index
                self.zones[nz].geometry.normal_coord, field_start = \
                    self._read_XYZblock(bytes_data, field_start, 3,
                                        self.zones[nz].iMax, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                # if file contains IBLANK (int) data, read that
                if self.has_iblank == True:
                    self.zones[nz].geometry.iblank, field_start = \
                        self._read_IBLANKblock(bytes_data, field_start,
                                               self.zones[nz].iMax,
                                               self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        else:
            # TODO: read non-constant geometry data
            print('Reading non-constant geometry data not implemented yet!')


    # ************************************************************************
    def _write_geometry_data(self, filename):
        """
        Writes PSU-WOPWOP geometry patch data to 'filename'
        Must be immediately preceded by a call to '_write_geometry_header(filename)'!
        """

        # open file in append mode - no need to adjust index
        with open(filename, 'ab') as f:

            # -----------------------------------------------------------
            # constant geometry
            if self.geometry_time_type == 'constant':

                # for each zone
                for nz in range(self.n_zones):

                    # write XYZ coords
                    for n in range(3):

                        # write order is Fortran (column-major)
                        for j in range(self.zones[nz].jMax):
                            for i in range(self.zones[nz].iMax):
                                _write_binary(f, self.zones[nz].geometry.XYZ_coord[n, i, j])

                    # write normal vector coords
                    for n in range(3):
                        for j in range(self.zones[nz].jMax):
                            for i in range(self.zones[nz].iMax):
                                _write_binary(f, self.zones[nz].geometry.normal_coord[n, i, j])

                    # write IBLANK data
                    if self.has_iblank == True:
                        for j in range(self.zones[nz].jMax):
                            for i in range(self.zones[nz].iMax):
                                _write_binary(f, self.zones[nz].geometry.iblank[i, j])

            # -----------------------------------------------------------
            else:
                # TODO: write non-constant geometry data
                print('Writing non-constant geometry data not implemented yet!')
            # -----------------------------------------------------------


    # **********************************************************************
    def _read_loading_header(self, filename):
        """
        Reads PSU-WOPWOP loading functional header info from 'filename'.
        """

        # do initial check for 'magic number' and endianness
        file_endianness = initial_check(filename)
        assert file_endianness == ENDIANNESS, \
            "File endianness does not match pywopwop default!"

        # open fileand read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # check version number - must be '1'and '0'
        self.version_number_major = _read_int(bytes_data, 4)
        self.version_number_minor = _read_int(bytes_data, 8)

        assert ((self.version_number_major == 1) and (self.version_number_minor == 0)), \
            'File version is not v1.0!'

        # read comments string (1024 bytes, starting at index 12)
        self.loading_comment = _read_string(bytes_data, 12, 1024)


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read format string (10 ints, 40 bytes, starting at index 1036),
        # assert loading file info vs. previously read from geometry file

        self.loading_format_string = []
        for n in range(10):
            self.loading_format_string.append(_read_int(bytes_data, 1036 + 4*n))

        # Verify file type description
        # --> loading_format_string[0] is '2' to indicate functional data file
        assert (self.loading_format_string[0] == 2), \
            "Format string does not start with '2', not a PSU-WOPWOP functional file!"

        assert (self.n_zones == self.loading_format_string[1]), \
            "Number of zones in loading file does not match geometry file!"

        assert (self.is_structured
                == _reverse_dict(structured_dict, self.loading_format_string[2])), \
            "Loading file 'is_structured' property does not match!"

        self.loading_time_type = \
            _reverse_dict(loading_time_dict, self.loading_format_string[3])

        assert (self.centered_type
                == _reverse_dict(centered_dict, self.loading_format_string[4])), \
            "Loading file 'centered_type' property does not match!"

        self.loading_data_type  = \
            _reverse_dict(loading_data_dict, self.loading_format_string[5])

        self.loading_ref_frame  = \
            _reverse_dict(ref_frame_dict, self.loading_format_string[6])

        assert (self.float_type
                == _reverse_dict(float_dict, self.loading_format_string[7])), \
            "Loading file 'float_type' property does not match!"

        # --> loading_format_string[8] is reserved for future use, and must be '0' in this version
        # --> loading_format_string[9] is reserved for future use, and must be '0' in this version


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # zone specification

        # read number of zones containing functional data (4 bytes)
        self.n_zones_with_loading_data = _read_int(bytes_data, 1076)

        # create list of zones containing functional data
        # -->>  negative numbers indicate zones for which WOPWOP should *NOT*
        #       calculate thickness noise - e.g. loading patch
        #
        # -->> PSU_WOPWOP zone list indices are one-based - subtract one to
        #       obtain zero-based (Python) indices
        self.zones_with_loading_data = []

        for z in range(self.n_zones_with_loading_data):
            self.zones_with_loading_data.append(_read_int(bytes_data, 1080 + z*4) - 1)


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read header info

        zone_info_start = 1080 + self.n_zones_with_loading_data*4

        # version 1.0
        if self.version_number_minor == 0:

            # --------------------------------------------------------------
            # structured, constant loading
            if self.is_structured == True:

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                if self.loading_time_type == 'constant':

                    for i, nz in enumerate(self.zones_with_loading_data):

                        # remove sign from zone index
                        nz = abs(nz)

                        # get handle to existing zone in list
                        zone = self.zones[nz]

                        # read loading zone name, assert iMax and jMax match
                        zone.loading_name = _read_string(bytes_data, zone_info_start + i*zone.header_length, 32)
                        iMax_fromfile = _read_int(bytes_data, zone_info_start + 32 + i*zone.header_length)
                        jMax_fromfile = _read_int(bytes_data, zone_info_start + 36 + i*zone.header_length)

                        assert ((zone.iMax == iMax_fromfile) and (zone.jMax == jMax_fromfile)), \
                            "(iMax, jMax) from loading file don't match existing values in PWWPatch instance!"

                        # set loading data flag
                        zone.has_loading_data = True

                        #self.zones.append(zone)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                else:
                    # TODO: implement non-constant loading
                    print('Reading non-constant functional zone info not implemented yet!')

            # ------------------------------------------------------------------
            else:
                # TODO: implement non-structured loading
                print('Reading non-structured functional zone info not implemented yet!')

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # TODO: implement header reader for version 1.1
        elif self.version_number2 == '1':
            print('Reading functional data file v1.1 header not implemented yet!')


    # ***********************************************************************
    def _write_loading_header(self, filename):
        """
        Writes PSU-WOPWOP functional file header to 'filename'
        """

        # build format string with most recent format values
        self._build_loading_format_string()

        with open(filename, 'wb') as f:

            # write 'magic number' 42 to first 4 bytes
            _write_binary(f, MAGICNUMBER)

            # write version number
            _write_binary(f, self.version_number_major)
            _write_binary(f, self.version_number_minor)

            # write comments string (1024 bytes)
            _write_string(f, self.loading_comment, 1024)

            # write format string (10 ints, 40 bytes)
            for n in range(10):
                _write_binary(f, self.loading_format_string[n])

            # --------------------------------------------------------------
            # write zone specification

            # write number of zones with data
            _write_binary(f, self.n_zones_with_loading_data)

            # write list of those zones - add one to create one-based
            # (PSU-WOPWOP) indices
            for z in self.zones_with_loading_data:
                _write_binary(f, (z + 1))

            # --------------------------------------------------------------
            # write zone info

            if self.is_structured == True:
                if self.loading_time_type == 'constant':

                    # for each zone containing data...
                    for nz in self.zones_with_loading_data:

                        # write name (32-byte string)
                        _write_string(f, self.zones[nz].loading_name, 32)

                        # write iMax and jMax (4 byte ints)
                        _write_binary(f, self.zones[nz].iMax)
                        _write_binary(f, self.zones[nz].jMax)

                else:
                    # TODO: implement non-constant functional data header
                    print('Writing non-constant functional data header not implemented yet!')

            # ------------------------------------------------------------
            else:
                # TODO: implement non-structured headers
                print('Writing non-structured functional data header not implemented yet!')


    # ***********************************************************************
    def _read_loading_data(self, filename):
        """
        Reads PSU-WOPWOP loading data from 'filename'.
        """

        # open file and read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # structured loading
        if self.is_structured == True:

            # ----------------------------------------------------------------
            if self.loading_time_type == 'constant':

                # start index for reading functional data
                # end of format string + zone specification + header
                field_start = (1076
                               + (1 + self.n_zones_with_loading_data)*4
                               + self.n_zones_with_loading_data*self.zones[0].header_length)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
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
                            self._read_XYZblock(bytes_data, field_start, 1,
                                                self.zones[nz].iMax,
                                                self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
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
                            self._read_XYZblock(bytes_data, field_start, 3,
                                                self.zones[nz].iMax,
                                                self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                elif self.data_type == 'flow_params':
                    print('Reading flow parameters not implemented yet!')


            # --------------------------------------------------------------
            else:
                # TODO: read non-constant functional data
                print('Reading non-constant functional data not implemented yet!')

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        else:
            # TODO: read non-structured functional data
            print('Reading non-structured functional data not implemented yet!')


    # ***********************************************************************
    def _write_loading_data(self, filename):
        """
        Writes PSU-WOPWOP functional data to 'filename'
        Must be immediately preceded by a call to '_write_header(filename)'!
        """

        # open file in append mode - no need to adjust index
        with open(filename, 'ab') as f:

            # -----------------------------------------------------------
            # structured loading
            if self.is_structured == True:

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                # constant geometry
                if self.loading_time_type == 'constant':

                    # ......................................................
                    # if data is surface pressure
                    if self.loading_data_type == 'surf_pressure':

                        # for each zone
                        for nz in self.zones_with_loading_data:

                            # remove negative sign
                            nz = abs(nz)

                            # write pressure data in Fortran (column-major) order
                            for j in range(self.zones[nz].jMax):
                                for i in range(self.zones[nz].iMax):
                                    _write_binary(f, self.zones[nz].loading.pressures[i, j])

                    # ......................................................
                    # if data is surface loading vectors
                    elif self.loading_data_type == 'surf_loading_vec':

                        # for each zone
                        for nz in self.zones_with_loading_data:

                            # remove negative sign
                            nz = abs(nz)

                            # write loading vectors in Fortran (column-major) order
                            for n in range(3):
                                for j in range(self.zones[nz].jMax):
                                    for i in range(self.zones[nz].iMax):
                                        _write_binary(f, self.zones[nz].loading.loading_vectors[n, i, j])

                    # ......................................................
                    elif self.loading_data_type == 'flow_params':
                        # TODO: write flow params data
                        print('Writing flow parameters not implemented yet!')

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                else:
                    # TODO: write non-constant loading data
                    print('Writing non-constant loading data not implemented yet!')

            # -----------------------------------------------------------
            else:
                # TODO: write non-structured loading data
                print('Writing non-structured loading data not implemented yet!')


    # **********************************************************************
    def _update_n_zones(self):
        self.n_zones = len(self.zones)
        self.n_zones_with_loading_data = len(self.zones_with_loading_data)


    # **********************************************************************
    def _build_loading_format_string(self):
        """
        Create list of data format values, later written to loading file as
        string of ints.
        """

        self._update_n_zones()

        self.loading_format_string = []
        self.loading_format_string.append(2)               # indicate functional data file
        self.loading_format_string.append(self.n_zones)
        self.loading_format_string.append(structured_dict[self.is_structured])
        self.loading_format_string.append(loading_time_dict[self.loading_time_type])
        self.loading_format_string.append(centered_dict[self.centered_type])
        self.loading_format_string.append(loading_data_dict[self.loading_data_type])
        self.loading_format_string.append(ref_frame_dict[self.loading_ref_frame])
        self.loading_format_string.append(float_dict[self.float_type])
        self.loading_format_string.append(RESERVED_DIGIT)   # reserved digit '0' inserted twice
        self.loading_format_string.append(RESERVED_DIGIT)


    def _build_geometry_format_string(self):
        """
        Create list of data format values, later written to geometry file as
        string of ints.
        """

        self._update_n_zones()

        self.geometry_format_string = []
        self.geometry_format_string.append(geom_dict[self.geometry_type])
        self.geometry_format_string.append(self.n_zones)
        self.geometry_format_string.append(structured_dict[self.is_structured])
        self.geometry_format_string.append(geometry_time_dict[self.geometry_time_type])
        self.geometry_format_string.append(centered_dict[self.centered_type])
        self.geometry_format_string.append(float_dict[self.float_type])
        self.geometry_format_string.append(iblank_dict[self.has_iblank])
        self.geometry_format_string.append(RESERVED_DIGIT)


# %% ##########################################################################
class Zone:
    """
    Parent class for zone data containing zone name, header length (fixed at
    32 bytes) and some calculation flags.
    """
    def __init__(self):
        self.name = ''
        self.number = 0                     # zone numbers use 1-based index

        self.calc_thickness_noise = True
        self.has_loading_data = False

        self.header_length = 32             # 'name' string has 32 bytes

        self.geometry = None
        self.loading = None


    def __str__(self):
        str1 = '\n\tGeometry name:          ' + self.geometry_name
        str2 = '\n\tLoading name:           ' + self.loading_name
        str3 = '\n\tGeometry info:          ' + self.geometry_info_str
        str4 = '\n\tCalc thickness noise:   ' + str(self.calc_thickness_noise)
        str5 = '\n\tHas loading data:       ' + str(self.has_loading_data)

        return str1+str2+str3+str4+str5


# ##########################################################################
class StructuredZone(Zone):
    """
    Parent class for structured zone, containing structured dimensions (2 ints)
    """

    def __init__(self):
        super().__init__()
        self.iMax = 0
        self.jMax = 0
        self.header_length += 2*VALUE_LENGTH

        self._update_geometry_info_str()


    def _update_geometry_info_str(self):
        str_iMax = '\n\t--> iMax:               ' + str(self.iMax)
        str_jMax = '\n\t--> jMax:               ' + str(self.jMax)
        self.geometry_info_str = str_iMax + str_jMax


    # **********************************************************************
    def add_StructuredConstantGeometry(self, XYZ_coord, normal_coord):
        """
        Adds structured, constant geometry data to current structured zone.

        Parameters
        ----------
        XYZ_coord : (3, iMax, jMax) array_like
            Array of mesh point coordinates to be added.

        normal_coord : (3, iMax, jMax) array_like
            Array of normal vector coordinates to be added.

        Returns
        -------
        None.
        """

        # updates iMax, jMax
        _, self.iMax, self.jMax = XYZ_coord.shape

        self.geometry = StructuredConstantGeometry(XYZ_coord, normal_coord)
        self._update_geometry_info_str()


    def add_StructuredPeriodicGeometry(self, XYZ_coord, normal_coord):
        # TODO: implement Structured Periodic Geometry
        print("Structured Periodic Geometry data not implemented yet!")


    def add_StructuredAperiodicGeometry(self, XYZ_coord, normal_coord):
        # TODO: implement Structured Aperiodic Geometry
        print("Structured Aperiodic Geometry data not implemented yet!")


    # **********************************************************************
    def add_StructuredConstantLoading(self, loading_data, loading_data_type):
        """
        Adds structured constant loading to current structured zone.

        Parameters
        ----------
        loading_data : (iMax, jMax) or (3, iMax, jMax) array_like
            The array of data to be added. Its shape is (iMax, jMax) for
            pressure data, and (3, iMax, jMax) for loading vector data.

        loading_data_type : {'surf_pressure', 'surf_loading_vec', 'flow_params'} string
            A string describing the type of loading data.

        Returns
        -------
        None.

        """

        # check 'loading_data_type' arg vs. loading data array shape
        if loading_data_type == 'surf_pressure':
            assert loading_data.shape == (self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for 'surf_pressure' (iMax, jMax)!"

        elif loading_data_type == 'surf_loading_vec':
            assert loading_data.shape == (3, self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for 'surf_loading_vec' (3, iMax, jMax)!"

        elif loading_data_type == 'flow_params':
            print("Cannot add 'flow_params' data to StructuredZone - not implemented yet!")

        self.loading = StructuredConstantLoading(loading_data, loading_data_type)


    # **********************************************************************
    def add_StructuredPeriodicLoading(self, loading_data, loading_data_type):
        # TODO: implement Structured Periodic Loading
        print("Structured Periodic Loading data not implemented yet!")


    # **********************************************************************
    def add_StructuredAperiodicLoading(self, loading_data, loading_data_type):
        # TODO: implement Structured Aperiodic Loading
        print("Structured Aperiodic Loading data not implemented yet!")


class StructuredConstantGeometry():
    """
    Class to store structured, constant geometry data - contains no time
    information.

    Parameters
    ----------
    XYZ_coord : (3, iMax, jMax) array_like
        Array of mesh point coordinates to be added.

    normal_coord : (3, iMax, jMax) array_like
        Array of normal vector coordinates to be added.

    Returns
    -------
    None.
    """

    def __init__(self, XYZ_coord, normal_coord):
        self.XYZ_coord = np.copy(XYZ_coord)
        self.normal_coord = np.copy(normal_coord)


class StructuredConstantLoading():
    """
    Class to store structured, constant pressure or loading vector data.

    Parameters
    ----------
    loading_data : (iMax, jMax) or (3, iMax, jMax) array_like
        The array of data to be added. Its shape is (iMax, jMax) for
        pressure data, and (3, iMax, jMax) for loading vector data.

    loading_data_type : {'surf_pressure', 'surf_loading_vec', 'flow_params'} string
        A string describing the type of loading data.

    Returns
    -------
    None.

    """

    def __init__(self, loading_data, loading_data_type):

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # if data is (iMax, jMax)-shaped array of surface pressures
        if loading_data_type == 'surf_pressure':

            # copy input data
            self.pressures = np.copy(loading_data)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # if data is (3, iMax, jMax)-shaped array of surface loading vectors
        elif loading_data_type == 'surf_loading_vec':

            # copy input data
            self.loading_vectors = np.copy(loading_data)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        elif loading_data_type == 'flow_params':
            # TODO: implement structured constant loading using flow params!
            print('Structured constant loading using flow params not implemented yet!')



# %% ##########################################################################
# class UnstructuredZone(Zone):
#     """
#     Parent class for unstructured zone, containing number of notes (int),
#     number of faces (int), and a connectivity list.

#     The connectivity list contains lists of integers defining each face: each
#     face is specified by one int indicating how many nodes it contains,
#     followed by that many integers defining the face. It is not used internally
#     by PSU-WOPWOP, but is used to output sigma surfaces.

#     Indices are one-based instead of zero-based, and are ordered clockwise
#     around the outward-facing normal for each face.
#     """

#     def __init__(self):
#         super().__init__()
#         self.nbNodes = 0
#         self.nbFaces = 0
#         self.connectivity = None
#         self.header_length += 2*VALUE_LENGTH

#         # TODO: define geometry_info_str!

# **********************************************************************
# class PeriodicZone(Zone):
#     """
#     Parent class for time-varying period zone, containing period (in seconds)
#     and number of time steps (integer)
#     """
#     def __init__(self):
#         super().__init__()
#         self.period = 0.
#         self.nt = 0
#         self.header_length += 2*VALUE_LENGTH


# class AperiodicZone(Zone):
#     """
#     Parent class for time-varying aperiod zone, containing number of time
#     steps (integer)
#     """
#     def __init__(self):
#         super().__init__()
#         self.nt = 0
#         self.header_length += VALUE_LENGTH


# # TODO: add period and nt info on initialization
# class StructuredPeriodicZone(StructuredZone, PeriodicZone):
#     """
#     Class for structured, periodic zone
#     """
#     def __init__(self):
#         super().__init__()


# # TODO: add nt info on initialization
# class StructuredAperiodicZone(StructuredZone, AperiodicZone):
#     """
#     Class for structured, aperiodic zone
#     """
#     def __init__(self):
#         super().__init__()


# %% ##########################################################################
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
        print('Magic number is correct - file {} is little endian\n'.format(filename))

    # if big endian, change flag and continue
    elif bytes_data == b'\x00\x00\x00*':
        endianness_flag = 'big'
        print('Magic number is correct - file {} is big endian\n'.format(filename))

    # if none, magic number is incorrect - probably not PSU-WOPWOP file!
    else:
        raise ValueError('Magic number is incorrect - file {} is probably not a PSU-WOPWOP patch file v1.0!'.format(filename))

    return endianness_flag


def _read_int(obj_name, start_index, n_bytes=VALUE_LENGTH,
              endianness_flag=ENDIANNESS):
    """
    Reads one integer value from an open file and returns the unpacked value.

    Parameters
    ----------
    obj_name : bytes_object
        Object containing binary data, such as an open file object.

    start_index : int
        Starting index of value to be read.

    n_bytes : {1, 2, 4, 8}
        Size of the integer to be read, in bytes. The default is the constant
        'VALUE_LENGTH' = 4.

    endianness_flag : {'little', 'big'}, optional
        String indicating the byte endinaness to be used. The default is
        the constant 'ENDIANNESS' = 'little'.

    Returns
    -------
    out : int
        Integer value unpacked from file data.

    """

    n_bytes_dict = {1:'b', 2:'h', 4:'i', 8:'q'}

    endianness_dict = {'little':'<', 'big':'>'}

    return struct.unpack(endianness_dict[endianness_flag] + n_bytes_dict[n_bytes],
                         obj_name[start_index:start_index + n_bytes])[0]


def _read_float(obj_name, start_index, n_bytes=VALUE_LENGTH,
                endianness_flag=ENDIANNESS):
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


def _write_binary(file, data, length=VALUE_LENGTH,
                  endianness_flag=ENDIANNESS, is_signed=IS_SIGNED):
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

