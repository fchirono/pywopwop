"""
pywopwop - https://github.com/fchirono/pywopwop

    Collection of convenience routines to parse and create PSU-WOPWOP input
    files version 1.0.

    --> PSU-WOPWOP Sigma surface file converter for visualization with
        Paraview.

Author:
    Fabio Casagrande Hirono
    Dec 2021
"""


import numpy as np

from readers_and_writers import read_block, write_block


def extract_var_names(nam_filename):
    """
    Parses a 'sigma.nam' file containing the variable names, and outputs a
    list of these names.
    """
    
    names = []
    with open(nam_filename, 'r') as file:
        for line in file:
            names.append(line.strip())
    
    return names