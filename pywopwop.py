# -*- coding: utf-8 -*-
"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

Author:
    Fabio Casagrande Hirono
    Dec 2021
"""


import numpy as np

from consts_and_dicts import MAGICNUMBER, ENDIANNESS, VALUE_LENGTH, \
    IS_SIGNED, RESERVED_DIGIT, reverse_dict, geom_dict, structured_dict, \
    loading_time_dict, geometry_time_dict, centered_dict, loading_data_dict, \
    ref_frame_dict, float_dict, iblank_dict

from readers_and_writers import initial_check, read_block, write_block, \
    read_IBLANKblock, read_int, read_float, write_binary, write_string, \
    read_string

from zones import Zone, StructuredZone, StructuredConstantGeometry, \
    StructuredConstantLoading

from sigma_processing import extract_sigma_var_names, process_sigma_fn_file, \
    process_sigma_geom_file, write_p3d_file, process_sigma_files


# %% #######################################################################
# PSU-WOPWOP main class PWWPatch
# ##########################################################################

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
    def _set_string(self, string, attr_name, length):
        """
        Writes input string to given attribute name, enforcing ASCII
        compatibility and string length.
        """

        # check string is ASCII compatible
        ascii_error = 'String is not ASCII compatible!'
        assert string[:length].isascii(), ascii_error

        # check string has maximum length, pad with spaces otherwise
        if len(string) < length:
            string += (length-len(string))*' '

        setattr(self, attr_name, string[:length])


    def set_units_string(self, string):
        self._set_string(string, 'units_string', 32)

    def set_geometry_comment(self, string):
        self._set_string(string, 'geometry_comment', 1024)

    def set_loading_comment(self, string):
        self._set_string(string, 'loading_comment', 1024)


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

            loading_data : (iMax, jMax) or (3, iMax, jMax) or (5, iMax, jMax) array_like
                Array of constant loading data to be added. Its shape is
                (iMax, jMax) for pressure data, (3, iMax, jMax) for loading
                vector data, and (5, iMax, jMax) for flow parameters
                (rho, rho*u, rho*v, rho*w, p').

        """

        # check PWWPatch instance is indeed structured
        assert (self.is_structured == True), \
            "Cannot add structured zone - PWWPatch instance is not structured!"

        # instantiate new zone
        zone = StructuredZone()

        zone.iMax, zone.jMax = XYZ_coord.shape[1:]

        # cap 'name' length to 32 bytes
        #zone.geometry_name = name[:32]
        #zone.loading_name = name[:32]
        zone.set_name(name)
        zone._set_string(name, 'geometry_name', 32)
        zone._set_string(name, 'loading_name', 32)


        # zone number will always be the current length of 'zones'
        zone.number = len(self.zones)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # adds geometry data to zone
        if self.geometry_time_type == 'constant':
            zone.add_StructuredConstantGeometry(XYZ_coord, normal_coord)

        elif self.geometry_time_type == 'periodic':
            # TODO: implement structured periodic geometry
            raise NotImplementedError("Can't add Periodic Geometry data to StructuredZone - not implemented yet!")

        elif self.geometry_time_type == 'aperiodic':
            # TODO: implement structured aperiodic geometry
            raise NotImplementedError("Can't add Aperiodic Geometry data to StructuredZone - not implemented yet!")

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
                    #print("Can't check loading data shape for flow parameter data - not implemented yet!")
                    assert loading_data.shape == (5, zone.iMax, zone.jMax), \
                        "'loading_data' does not match expected shape for 'flow_params' (5, iMax, jMax)!"

                zone.add_StructuredConstantLoading(loading_data, self.loading_data_type)

            # ----------------------------------------------------------------
            elif self.loading_time_type == 'periodic':

                # TODO: implement structured periodic loading
                raise NotImplementedError("Can't add Structured Periodic Loading data to StructuredZone - not implemented yet!")

            # ----------------------------------------------------------------
            elif self.loading_time_type == 'aperiodic':

                # TODO: implement structured aperiodic loading
                raise NotImplementedError("Can't add Structured Aperiodic Loading data to StructuredZone - not implemented yet!")

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
        self.version_number_major = read_int(bytes_data, 4)
        self.version_number_minor = read_int(bytes_data, 8)

        assert ((self.version_number_major == 1)
                and (self.version_number_minor == 0)), \
            'File version is not v1.0!'

        # read units string (32 chars, starting at index 12)
        self.units_string = read_string(bytes_data, 12, 32)

        # read comments string (1024 bytes, starting at index 44)
        self.geometry_comment = read_string(bytes_data, 44, 1024)

        # read format string (8 ints, 32 bytes, starting at index 1068)
        self.geometry_format_string = []
        for n in range(8):
            self.geometry_format_string.append(read_int(bytes_data, 1068 + 4*n))

        # Populate file type description
        self.geometry_type      = reverse_dict(geom_dict, self.geometry_format_string[0])
        self.n_zones            = self.geometry_format_string[1]
        self.is_structured      = reverse_dict(structured_dict, self.geometry_format_string[2])
        self.geometry_time_type = reverse_dict(geometry_time_dict, self.geometry_format_string[3])
        self.centered_type      = reverse_dict(centered_dict, self.geometry_format_string[4])
        self.float_type         = reverse_dict(float_dict, self.geometry_format_string[5])
        self.has_iblank         = reverse_dict(iblank_dict, self.geometry_format_string[6])


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read zone info (starting at index 1100)

        if self.is_structured == True:

            if self.geometry_time_type == 'constant':

                for nz in range(self.n_zones):
                    # instantiate zone and read info from file
                    zone = StructuredZone()

                    zone.number = len(self.zones)

                    # reads geometry zone name
                    name = read_string(bytes_data, 1100 + nz*zone.header_length, 32)
                    zone.set_name(name)
                    zone.geometry_name = name
                    zone.loading_name = ''

                    zone.iMax = read_int(bytes_data, 1100 + 32 + nz*zone.header_length)
                    zone.jMax = read_int(bytes_data, 1100 + 36 + nz*zone.header_length)

                    zone._update_geometry_info_str()
                    self.zones.append(zone)

            else:
                # TODO: implement non-constant geometries
                raise NotImplementedError("Can't read non-constant geometry zone info - not implemented yet!")

        else:
            # TODO: implement non-structured headers
            raise NotImplementedError("Can't read non-structured geometry zone info - not implemented yet!")

        assert (self.n_zones == len(self.zones)), \
            "Number of zones in format string doesn't match file data!"


    # **********************************************************************
    def _write_geometry_header(self, filename):
        """
        Writes PSU-WOPWOP geometry header to 'filename'.
        """

        # build format string with most recent format values
        self._build_geometry_format_string()

        # open or create new file
        with open(filename, 'wb+') as f:

            # write 'magic number' 42 to first 4 bytes
            write_binary(f, MAGICNUMBER)

            # write version number
            write_binary(f, self.version_number_major)
            write_binary(f, self.version_number_minor)

            # write units string (32 chars)
            write_string(f, self.units_string, 32)

            # write comments string (1024 bytes)
            write_string(f, self.geometry_comment, 1024)

            # write format string (8 ints, 32 bytes)
            for n in range(8):
                write_binary(f, self.geometry_format_string[n])

            # --------------------------------------------------------------
            # write zone info
            if self.is_structured == True:

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                if self.geometry_time_type == 'constant':

                    # for each zone...
                    for nz in range(self.n_zones):

                        # write name (32-byte string)
                        write_string(f, self.zones[nz].geometry_name, 32)

                        # write iMax and jMax (4 byte ints)
                        write_binary(f, self.zones[nz].iMax)
                        write_binary(f, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                else:
                    # TODO: implement non-constant geometries
                    raise NotImplementedError("Can't write non-constant geometry zone info - not implemented yet!")

            # --------------------------------------------------------------
            else:
                # TODO: implement non-structured headers
                raise NotImplementedError("Can't write non-structured geometry zone info - not implemented yet!")


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
                    read_block(bytes_data, field_start, 3,
                               self.zones[nz].iMax, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                # read normal vector coords and next index
                self.zones[nz].geometry.normal_coord, field_start = \
                    read_block(bytes_data, field_start, 3,
                               self.zones[nz].iMax, self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                # if file contains IBLANK (int) data, read that
                if self.has_iblank == True:
                    self.zones[nz].geometry.iblank, field_start = \
                        read_IBLANKblock(bytes_data, field_start,
                                         self.zones[nz].iMax,
                                         self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        else:
            # TODO: read non-constant geometry data
            raise NotImplementedError("Can't read non-constant geometry data - not implemented yet!")


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

                    write_block(f, self.zones[nz].geometry.XYZ_coord)
                    write_block(f, self.zones[nz].geometry.normal_coord)

                    if self.has_iblank == True:
                        write_block(f, self.zones[nz].geometry.iblank)

            # -----------------------------------------------------------
            else:
                # TODO: write non-constant geometry data
                raise NotImplementedError("Can't write non-constant geometry data - not implemented yet!")
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

        # open file and read binary content
        with open(filename, 'rb') as f:
            bytes_data = f.read()

        # check version number - must be '1'and '0'
        self.version_number_major = read_int(bytes_data, 4)
        self.version_number_minor = read_int(bytes_data, 8)

        assert ((self.version_number_major == 1) and (self.version_number_minor == 0)), \
            'File version is not v1.0!'

        # read comments string (1024 bytes, starting at index 12)
        self.loading_comment = read_string(bytes_data, 12, 1024)


        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # read format string (10 ints, 40 bytes, starting at index 1036),
        # assert loading file info vs. previously read from geometry file

        self.loading_format_string = []
        for n in range(10):
            self.loading_format_string.append(read_int(bytes_data, 1036 + 4*n))

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

                        # read loading zone name, assert iMax and jMax match
                        zone.loading_name = read_string(bytes_data, zone_info_start + i*zone.header_length, 32)
                        iMax_fromfile = read_int(bytes_data, zone_info_start + 32 + i*zone.header_length)
                        jMax_fromfile = read_int(bytes_data, zone_info_start + 36 + i*zone.header_length)

                        assert ((zone.iMax == iMax_fromfile) and (zone.jMax == jMax_fromfile)), \
                            "(iMax, jMax) from loading file don't match existing values in PWWPatch instance!"

                        # set loading data flag
                        zone.has_loading_data = True

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                else:
                    # TODO: implement non-constant loading
                    raise NotImplementedError("Can't read non-constant loading data - not implemented yet!")

            # ------------------------------------------------------------------
            else:
                # TODO: implement non-structured loading
                raise NotImplementedError("Can't read non-structured loading data - not implemented yet!")

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # TODO: implement header reader for version 1.1
        elif self.version_number2 == '1':
            raise NotImplementedError("Can't read loading data file v1.1 header - not implemented yet!")



    # ***********************************************************************
    def _write_loading_header(self, filename):
        """
        Writes PSU-WOPWOP functional file header to 'filename'
        """

        # build format string with most recent format values
        self._build_loading_format_string()

        with open(filename, 'wb') as f:

            # write 'magic number' 42 to first 4 bytes
            write_binary(f, MAGICNUMBER)

            # write version number
            write_binary(f, self.version_number_major)
            write_binary(f, self.version_number_minor)

            # write comments string (1024 bytes)
            write_string(f, self.loading_comment, 1024)

            # write format string (10 ints, 40 bytes)
            for n in range(10):
                write_binary(f, self.loading_format_string[n])

            # --------------------------------------------------------------
            # write zone specification

            # write number of zones with data
            write_binary(f, self.n_zones_with_loading_data)

            # write list of those zones - add one to create one-based
            # (PSU-WOPWOP) indices
            for z in self.zones_with_loading_data:
                write_binary(f, (z + 1))

            # --------------------------------------------------------------
            # write zone info

            if self.is_structured == True:
                if self.loading_time_type == 'constant':

                    # for each zone containing data...
                    for nz in self.zones_with_loading_data:

                        # write name (32-byte string)
                        write_string(f, self.zones[nz].loading_name, 32)

                        # write iMax and jMax (4 byte ints)
                        write_binary(f, self.zones[nz].iMax)
                        write_binary(f, self.zones[nz].jMax)

                else:
                    # TODO: implement non-constant functional data header
                    raise NotImplementedError("Can't write non-constant loading data header - not implemented yet!")

            # ------------------------------------------------------------
            else:
                # TODO: implement non-structured headers
                raise NotImplementedError("Can't write non-structured loading data header - not implemented yet!")

            # ------------------------------------------------------------

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
                            read_block(bytes_data, field_start, 1,
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
                            read_block(bytes_data, field_start, 3,
                                       self.zones[nz].iMax,
                                       self.zones[nz].jMax)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                elif self.loading_data_type == 'flow_params':
                    # TODO: read flow_params loading data
                    #raise NotImplementedError("Can't read flow_params loading data - not implemented yet!")

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

            # --------------------------------------------------------------
            else:
                # TODO: read non-constant loading data
                raise NotImplementedError("Can't read non-constant loading data - not implemented yet!")

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        else:
            # TODO: read non-structured loading data
            raise NotImplementedError("Can't read non-structured loading data - not implemented yet!")


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
                    if self.loading_data_type == 'surf_pressure':

                        # for each zone...
                        for nz in self.zones_with_loading_data:
                            nz = abs(nz)
                            write_block(f, self.zones[nz].loading.pressures)

                    # ......................................................
                    elif self.loading_data_type == 'surf_loading_vec':

                        # for each zone...
                        for nz in self.zones_with_loading_data:
                            nz = abs(nz)
                            write_block(f, self.zones[nz].loading.loading_vectors)

                    # ......................................................
                    elif self.loading_data_type == 'flow_params':
                        # TODO: write flow params data
                        #raise NotImplementedError("Can't write flow_params loading data - not implemented yet!")

                        # for each zone...
                        for nz in self.zones_with_loading_data:
                            nz = abs(nz)
                            write_block(f, self.zones[nz].loading.flow_params)

                # -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                else:
                    # TODO: write non-constant loading data
                    raise NotImplementedError("Can't write non-constant loading data - not implemented yet!")

            # -----------------------------------------------------------
            else:
                # TODO: write non-structured loading data
                raise NotImplementedError("Can't write non-structured loading data - not implemented yet!")


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
        self.loading_format_string.append(2)    # indicate functional data file
        self.loading_format_string.append(self.n_zones)
        self.loading_format_string.append(structured_dict[self.is_structured])
        self.loading_format_string.append(loading_time_dict[self.loading_time_type])
        self.loading_format_string.append(centered_dict[self.centered_type])
        self.loading_format_string.append(loading_data_dict[self.loading_data_type])
        self.loading_format_string.append(ref_frame_dict[self.loading_ref_frame])
        self.loading_format_string.append(float_dict[self.float_type])

        # reserved digit '0' inserted twice
        self.loading_format_string.append(RESERVED_DIGIT)
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


# %% #######################################################################
# Functions to compare two PWWPatch instances for identical content
# ##########################################################################

def _list_compare_attrs(obj1, obj2, attrs_to_ignore=[], level=0):
    """
    Lists and compares attributes from two objects, check them for identical
    content
    """

    is_equal = True

    # list all attributes in obj
    list_attrs = list(obj1.__dict__.keys())
    # remove attrs in 'attrs_to_ignore'
    for attr in attrs_to_ignore:
        list_attrs.remove(attr)

    for attr in list_attrs:
        attr1 = getattr(obj1, attr)
        attr2 = getattr(obj2, attr)

        if type(attr1) is np.ndarray:
            if not np.array_equal(attr1, attr2):
                print(level*'\t' + '{} : different'.format(attr))
                is_equal = False
        else:
            if attr1 != attr2:
                print(level*'\t' + '{} : different'.format(attr))
                is_equal = False

    return is_equal


def compare_pwwpatches(pwwpatch1, pwwpatch2):
    """
    Compares the contents of two PWWPatches, informs whether contents are
    identical or not.
    """

    is_equal = True

    # iterate and compare over both patches
    is_equal = _list_compare_attrs(pwwpatch1, pwwpatch2, ['zones',], 0)

    # if no differences were detected so far, check zones' contents as well
    if is_equal:
        print("Patches' attributes equal so far; comparing zones now...\t")

        for z in range(pwwpatch1.n_zones):
            print('Zone {}:'.format(z))
            zones_are_equal = _list_compare_attrs(pwwpatch1.zones[z],
                                                  pwwpatch2.zones[z],
                                                  ['geometry', 'loading'], 1)

            # compare geometry and loading contents
            geom_are_equal = _list_compare_attrs(pwwpatch1.zones[z].geometry,
                                                 pwwpatch2.zones[z].geometry,
                                                 [], 2)

            loading_is_equal = _list_compare_attrs(pwwpatch1.zones[z].loading,
                                                   pwwpatch2.zones[z].loading,
                                                   [], 2)

    is_equal = is_equal & zones_are_equal & geom_are_equal & loading_is_equal

    if is_equal:
        print('File contents are identical!')
    else:
        print('File contents are NOT identical!')

    return
