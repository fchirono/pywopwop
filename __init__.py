# -*- coding: utf-8 -*-
"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    The main class 'PWWPatch' is defined in '__init__.py', and imports methods
    from other files in the package.

Author:
    Fabio Casagrande Hirono
    Nov 2022
"""


import numpy as np


from ._zones import Zone, StructuredZone, StructuredConstantGeometry, \
    StructuredConstantLoading, StructuredAperiodicLoading

from ._geometry_readers_writers import _read_geometry_header, \
    _read_geometry_data, _write_geometry_header, _write_geometry_data

from ._loading_readers_writers import _read_loading_header, _read_loading_data,\
    _write_loading_header, _write_loading_data

from ._binary_readers_writers import initial_check, read_block, write_block, \
    read_IBLANKblock, read_int, read_float, write_binary, write_string, \
    read_string

from ._consts_and_dicts import MAGICNUMBER, ENDIANNESS, VALUE_LENGTH, \
    IS_SIGNED, RESERVED_DIGIT, reverse_dict, geom_dict, structured_dict, \
    loading_time_dict, geometry_time_dict, centered_dict, loading_data_dict, \
    ref_frame_dict, float_dict, iblank_dict

from ._sigma_processing import extract_sigma_var_names, process_sigma_fn_file, \
    process_sigma_geom_file, write_p3d_file, process_sigma_files



# ##########################################################################
# %%PSU-WOPWOP main class PWWPatch
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
    # file readers
    def read_geometry_file(self, geometry_filename):
        _read_geometry_header(self, geometry_filename)
        _read_geometry_data(self, geometry_filename)


    def read_loading_file(self, loading_filename):
        _read_loading_header(self, loading_filename)
        _read_loading_data(self, loading_filename)


    # file writers
    def write_geometry_file(self, geometry_filename):
        _write_geometry_header(self, geometry_filename)
        _write_geometry_data(self, geometry_filename)


    def write_loading_file(self, loading_filename):
        _write_loading_header(self, loading_filename)
        _write_loading_data(self, loading_filename)


    # **********************************************************************
    def add_StructuredZone(self, name, XYZ_coord, normal_coord,
                           calc_thickness_noise=True, loading_data=None,
                           time_steps=None):
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

        time_steps : (Nt,)-shaped array_like, optional
            Array containing the time steps for periodic or aperiodic data.
            If aperiodic, the values must be a time in seconds. If periodic,
            the values can be a time in seconds, or an angle in degrees, or
            another convenient measure.

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

        If  'geometry_time_type' is 'aperiodic', the geometry arrays are:

            XYZ_coord : (Nt, 3, iMax, jMax) array_like
                Array of mesh point coordinates to be added per timestep.

            normal_coord : (Nt, 3, iMax, jMax) array_like
                Array of normal vector coordinates to be added per timestep.


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

        # # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # # read mesh dimensions - already done when reading
        # if self.geometry_time_type == 'constant':
        #     zone.iMax, zone.jMax = XYZ_coord.shape[1:]

        # elif self.geometry_time_type == 'aperiodic':
        #     zone.Nt, _, zone.iMax, zone.jMax = XYZ_coord.shape

        #     # assert 'Nt' in geometry array is the same as in 'time_vec'
        #     assert (zone.Nt == time_steps.shape[0]), \
        #         "Current zone does not contain the same number of time steps as time vector!"

        # # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

        zone._set_string(name, 'name', 32)
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
            zone.add_StructuredAperiodicGeometry(XYZ_coord, normal_coord)

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


# ##########################################################################
# %% Functions to compare two PWWPatch instances for identical content
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

