"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    --> PSU-WOPWOP file readers and writers

Author:
    Fabio Casagrande Hirono
    Dec 2021
"""


import numpy as np
import struct

from consts_and_dicts import ENDIANNESS, VALUE_LENGTH, IS_SIGNED


# %% #######################################################################
# PSU-WOPWOP initial file check for 'magic number' and endianness
# ##########################################################################

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

    # if data is 42 in little endian, continue
    if bytes_data == b'*\x00\x00\x00':
        print('Magic number is correct - file {} is little endian\n'.format(filename))

    # if data is 42 in little endian, change flag and continue
    elif bytes_data == b'\x00\x00\x00*':
        endianness_flag = 'big'
        print('Magic number is correct - file {} is big endian\n'.format(filename))

    # if magic number is incorrect, it's probably not PSU-WOPWOP file!
    else:
        raise ValueError('Magic number is incorrect - file {} is probably not a PSU-WOPWOP patch file v1.0!'.format(filename))

    return endianness_flag


# %% #######################################################################
# PSU-WOPWOP Plot3D-like block readers
# ##########################################################################

def read_block(bytes_data, start_index, num_dims, iMax, jMax):
    """
    Reads a block of data in PLOT3D-like format from a binary file.

    The block of data is a (num_dims, iMax, jMax)-shaped array of float32
    values. It can be a block of 3D coordinates (x,y,z) if 'num_dims=3', 2D
    coordinates (x,y) if 'num_dims=2', or one-dimensional

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
    block_data : (num_dims, iMax, jMax) array_like
        Numpy array containing the block data using float32 numbers.

    next_index : int
        Index for the next data field in the binary file.


    Notes
    -----
    The data in the file is assumed to be organized in Fortran order:

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
        block_data = np.zeros((iMax, jMax), dtype=np.float32)

        for j in range(jMax):
            for i in range(iMax):

                # fields are (i,j) order, hence j*iMax + i
                current_index = (start_index
                                 + ((j*iMax + i)*VALUE_LENGTH))

                block_data[i, j] = read_float(bytes_data, current_index)

    # read XYZ or XY - i.e. (num_dims, iMax, jMax)-shaped
    else:
        block_data = np.zeros((num_dims, iMax, jMax), dtype=np.float32)

        for n in range(num_dims):
            for j in range(jMax):
                for i in range(iMax):

                    # fields are (i,j) order, hence j*iMax + i
                    current_index = (start_index
                                     + ((n*iMax*jMax + j*iMax + i)
                                        * VALUE_LENGTH))

                    block_data[n, i, j] = read_float(bytes_data, current_index)

    # increase start index by ('value_length' bytes * num_dims coords * iMax * jMax)
    next_index = start_index + VALUE_LENGTH*(num_dims*iMax*jMax)

    return block_data, next_index


def read_IBLANKblock(bytes_data, start_index, iMax, jMax):
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
                              + (j*iMax + i)*VALUE_LENGTH)

            IBLANK_data[i, j] = read_int(bytes_data, current_index)

    # increase start index by ('value_length' bytes * iMax * jMax)
    next_index = start_index + VALUE_LENGTH*(iMax*jMax)

    return IBLANK_data, next_index


# %% #######################################################################
# PSU-WOPWOP string reader / writer
# ##########################################################################

def write_string(file, string, max_length):
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


def read_string(obj_name, start_index, len_string):
    """
    Reads strings of arbitrary length from binary file object.
    """
    mystring = ''
    for i in range(len_string):
        mystring += chr(read_int(obj_name, start_index + i, 1))

    return mystring


# %% #######################################################################
# PSU-WOPWOP int and float reader / writer
# ##########################################################################

def read_int(obj_name, start_index, n_bytes=VALUE_LENGTH,
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


def read_float(obj_name, start_index, n_bytes=VALUE_LENGTH,
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


def write_binary(file, data, length=VALUE_LENGTH,
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
