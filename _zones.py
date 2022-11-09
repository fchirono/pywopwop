"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    --> PSU-WOPWOP zones class definition

Author:
    Fabio Casagrande Hirono
    Nov 2022
"""


import numpy as np

from ._consts_and_dicts import VALUE_LENGTH


# ##########################################################################
# %% PSU-WOPWOP parent class Zone
# ##########################################################################

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

        self.geometry_header_length = 32    # 'name' string has 32 bytes
        self.loading_header_length = 32

        self.geometry = None
        self.loading = None


    def __str__(self):
        str1 = '\n\tGeometry name:          ' + self.geometry_name
        str2 = '\n\tLoading name:           ' + self.loading_name
        str3 = '\n\tGeometry info:          ' + self.geometry_info_str
        str4 = '\n\tCalc thickness noise:   ' + str(self.calc_thickness_noise)
        str5 = '\n\tHas loading data:       ' + str(self.has_loading_data)

        return str1+str2+str3+str4+str5


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


# ##########################################################################
# %% PSU-WOPWOP classes for structured zone
# ##########################################################################

class StructuredZone(Zone):
    """
    Parent class for structured zone, containing structured dimensions (2 ints)
    """

    def __init__(self):
        super().__init__()
        self.iMax = 0
        self.jMax = 0
        self.geometry_header_length += 2*VALUE_LENGTH
        self.loading_header_length += 2*VALUE_LENGTH

        self._update_geometry_info_str()


    def _update_geometry_info_str(self):
        str_iMax = '\n\t--> iMax:               ' + str(self.iMax)
        str_jMax = '\n\t--> jMax:               ' + str(self.jMax)

        self.geometry_info_str = str_iMax + str_jMax

        # if zone has time information, display that too
        if hasattr(self, 'Nt'):
            str_Nt = '\n\t--> Nt:                 ' + str(self.Nt)
            self.geometry_info_str += str_Nt


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


    # **********************************************************************
    def add_StructuredPeriodicGeometry(self, XYZ_coord, normal_coord):
        # TODO: implement Structured Periodic Geometry
        raise NotImplementedError("Can't add Structured Periodic Geometry data - not implemented yet!")


    # **********************************************************************
    def add_StructuredAperiodicGeometry(self, XYZ_coord, normal_coord):
        """
        Adds structured, aperiodic geometry data to current structured zone.

        Parameters
        ----------
        XYZ_coord : (Nt, 3, iMax, jMax) array_like
            Array of mesh point coordinates to be added per timestep.

        normal_coord : (Nt, 3, iMax, jMax) array_like
            Array of normal vector coordinates to be added per timestep.

        Returns
        -------
        None.
        """

        # updates iMax, jMax
        self.Nt, _, self.iMax, self.jMax = XYZ_coord.shape

        # increase geometry_header_length (must contain 'Nt' as well)
        self.geometry_header_length += VALUE_LENGTH

        self.geometry = StructuredAperiodicGeometry(XYZ_coord, normal_coord)
        self._update_geometry_info_str()


    # **********************************************************************
    def add_StructuredConstantLoading(self, loading_data, loading_data_type):
        """
        Adds structured constant loading to current structured zone.

        Parameters
        ----------
        loading_data : (iMax, jMax) or (3, iMax, jMax) or (5, iMax, jMax) array_like
            The array of data to be added. Its shape is (iMax, jMax) for
            pressure data, (3, iMax, jMax) for loading vector data, and
            (5, iMax, jMax) for flow parameters (rho, rho*u, rho*v, rho*w, p').

        loading_data_type : {'surf_pressure', 'surf_loading_vec', 'flow_params'} string
            A string describing the type of loading data.

        Returns
        -------
        None.

        """

        # check 'loading_data_type' arg vs. loading data array shape
        if loading_data_type == 'surf_pressure':
            assert loading_data.shape == (self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for constant 'surf_pressure' (iMax, jMax)!"

        elif loading_data_type == 'surf_loading_vec':
            assert loading_data.shape == (3, self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for constant 'surf_loading_vec' (3, iMax, jMax)!"

        elif loading_data_type == 'flow_params':
            assert loading_data.shape == (5, self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for constant 'flow_params' (5, iMax, jMax)!"

        self.loading = StructuredConstantLoading(loading_data, loading_data_type)


    # **********************************************************************
    def add_StructuredPeriodicLoading(self, loading_data, loading_data_type):
        # TODO: implement Structured Periodic Loading
        raise NotImplementedError("Can't add Structured Periodic Loading data - not implemented yet!")


    # **********************************************************************
    def add_StructuredAperiodicLoading(self, loading_data, loading_data_type):

        # check 'loading_data_type' arg vs. loading data array shape
        if loading_data_type == 'surf_pressure':
            assert loading_data.shape == (self.Nt, self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for aperiodic 'surf_pressure' (Nt, iMax, jMax)!"

        elif loading_data_type == 'surf_loading_vec':
            assert loading_data.shape == (self.Nt, 3, self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for aperiodic 'surf_loading_vec' (Nt, 3, iMax, jMax)!"

        elif loading_data_type == 'flow_params':
            assert loading_data.shape == (self.Nt, 5, self.iMax, self.jMax), \
                "'loading_data' does not match expected shape for aperiodic 'flow_params' (Nt, 5, iMax, jMax)!"

        # increase loading_header_length (must contain 'Nt' as well)
        self.loading_header_length += VALUE_LENGTH

        self.loading = StructuredAperiodicLoading(loading_data, loading_data_type)



# ##########################################################################
# %% PSU-WOPWOP classes for structured geometries
# ##########################################################################

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

        assert XYZ_coord.ndim == 3, \
            "'XYZ_coord' dimensions do not match for Structured Constant Geometry data!"

        assert normal_coord.ndim == 3, \
            "'normal_coord' dimensions do not match for Structured Constant Geometry data!"

        self.XYZ_coord = np.copy(XYZ_coord)
        self.normal_coord = np.copy(normal_coord)


class StructuredAperiodicGeometry():
    """
    Class to store structured, aperiodic geometry data, containing the number
    of timesteps.

    Parameters
    ----------
    XYZ_coord : (Nt, 3, iMax, jMax) array_like
        Array of mesh point coordinates to be added at each timestep.

    normal_coord : (Nt, 3, iMax, jMax) array_like
        Array of normal vector coordinates to be added at each timestep.

    Returns
    -------
    None.
    """

    def __init__(self, XYZ_coord, normal_coord, Nt):

        assert XYZ_coord.ndim == 4, \
            "'XYZ_coord' dimensions do not match for Structured Aperiodic Geometry data!"

        assert normal_coord.ndim == 4, \
            "'normal_coord' dimensions do not match for Structured Aperiodic Geometry data!"


        self.XYZ_coord = np.copy(XYZ_coord)
        self.normal_coord = np.copy(normal_coord)

        # vector to store time values for each timestep
        self.t = np.zeros(Nt, dtype=np.float32)


# ##########################################################################
# %% PSU-WOPWOP classes for structured loading
# ##########################################################################

class StructuredConstantLoading():
    """
    Class to store structured, constant pressure, loading vector, or flow
    parameter data.

    Parameters
    ----------
    loading_data : (iMax, jMax) or (3, iMax, jMax) or (5, iMax, jMax) array_like
        The array of data to be added. Its shape is (iMax, jMax) for
        pressure data, (3, iMax, jMax) for loading vector data, and
        (5, iMax, jMax) for flow parameters (rho, rho*u, rho*v, rho*w, p').

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

            # assert data ndims
            assert loading_data.ndim == 2,\
                "'loading_data' dimensions do not match structured constant surface pressure!"

            # copy input data
            self.pressures = np.copy(loading_data)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # if data is (3, iMax, jMax)-shaped array of surface loading vectors
        elif loading_data_type == 'surf_loading_vec':

            # assert data ndims
            assert loading_data.ndim == 3,\
                "'loading_data' dimensions do not match structured constant surface loading vectors!"

            # copy input data
            self.loading_vectors = np.copy(loading_data)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # if data is (5, iMax, jMax)-shaped array of flow parameters
        elif loading_data_type == 'flow_params':

            # assert data ndims
            assert loading_data.ndim == 3,\
                "'loading_data' dimensions do not match structured constant flow parameters!"

            # copy input data
            self.flow_params = np.copy(loading_data)


class StructuredAperiodicLoading():
    """
    Class to store structured, aperiodic pressure, loading vector, or flow
    parameter data.

    Parameters
    ----------
    loading_data : (Nt, iMax, jMax) or (Nt, 3, iMax, jMax) or (Nt, 5, iMax, jMax) array_like
        The array of data to be added. Its shape is (Nt, iMax, jMax) for
        pressure data, (Nt, 3, iMax, jMax) for loading vector data, and
        (Nt, 5, iMax, jMax) for flow parameters (rho, rho*u, rho*v, rho*w, p').

    loading_data_type : {'surf_pressure', 'surf_loading_vec', 'flow_params'} string
        A string describing the type of loading data.

    Returns
    -------
    None.

    """

    def __init__(self, loading_data, loading_data_type):

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # if data is (Nt, iMax, jMax)-shaped array of surface pressures
        if loading_data_type == 'surf_pressure':

            # assert data ndims
            assert loading_data.ndim == 3,\
                "'loading_data' dimensions do not match structured aperiodic surface pressure!"

            # copy input data
            self.pressures = np.copy(loading_data)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # if data is (Nt, 3, iMax, jMax)-shaped array of surface loading vectors
        elif loading_data_type == 'surf_loading_vec':

            # assert data ndims
            assert loading_data.ndim == 4,\
                "'loading_data' dimensions do not match structured aperiodic surface loading vectors!"

            # copy input data
            self.loading_vectors = np.copy(loading_data)

        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        # if data is (Nt, 5, iMax, jMax)-shaped array of flow parameters
        elif loading_data_type == 'flow_params':

            # assert data ndims
            assert loading_data.ndim == 4,\
                "'loading_data' dimensions do not match structured aperiodic flow parameters!"

            # copy input data
            self.flow_params = np.copy(loading_data)

