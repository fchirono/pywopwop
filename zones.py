"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    --> PSU-WOPWOP zones class definition

Author:
    Fabio Casagrande Hirono
    Dec 2021
"""


import numpy as np

from consts_and_dicts import VALUE_LENGTH


# %% #######################################################################
# PSU-WOPWOP parent class Zone
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


# %% #######################################################################
# PSU-WOPWOP classes for structured data
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
        raise NotImplementedError("Structured Periodic Geometry data not implemented yet!")


    def add_StructuredAperiodicGeometry(self, XYZ_coord, normal_coord):
        # TODO: implement Structured Aperiodic Geometry
        raise NotImplementedError("Structured Aperiodic Geometry data not implemented yet!")


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
        raise NotImplementedError("Structured Periodic Loading data not implemented yet!")


    # **********************************************************************
    def add_StructuredAperiodicLoading(self, loading_data, loading_data_type):
        # TODO: implement Structured Aperiodic Loading
        raise NotImplementedError("Structured Aperiodic Loading data not implemented yet!")


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
            raise NotImplementedError('Structured Constant loading using flow params not implemented yet!')



# %% #######################################################################
# PSU-WOPWOP classes for unstructured data - UNDER DEVELOPMENT
# ##########################################################################

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
