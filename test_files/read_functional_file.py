# -*- coding: utf-8 -*-
"""
Reads PSU-WOPWOP binary functional (.dat) files using pywopwop

File contents can be compared in Cygwin using:
    cmp ./path-to/file1.dat ./path-to/file2.dat

Author:
    Fabio Casagrande Hirono
    Sep 2021
"""


import pywopwop as PWW
import matplotlib.pyplot as plt


# %% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# Opens functional files for inspection

# case 1 - constant geometry, constant compact loading
path = '../PSU-WOPWOP_v3.4.4/case1/'
filename_loading = 'gyrodyneLoading.dat'


# # case 5 - constant geometry, aperiodic loading
# path = '../PSU-WOPWOP_v3.4.4/case5/'
# filename_loading = 'AperLoadingShort.dat'

# %%
# create instance of  loading functional file
case_loadingfile = PWW.LoadingFunctionalFile(path + filename_loading)


# writes functional data to different functional file
filename2 = 'MyNewFile_loading.dat'
case_loadingfile.write_functional_file(filename2)
