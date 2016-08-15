"""
Functions for analyzing ion coordination in PDB structures
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import os.path

def en(protein, ion, maxdistance=20, oxynotprotein=True, periodic=True):
    """Gives the distances of oxygen atoms from an ion.
    :Arguments:
        *protein*
            protein Universe
        *ion*
            ion Atom
        *maxdistance*
            maximum distance of interest from the ion; default = 20
        *oxynotprotein*
            boolean value of whether to include oxygens not in the protein; default = True

    :Returns:
        *df*
            `pandas.DataFrame` containing resids, resnames, and atom names
            for each oxygen in the protein file
    """
    columns = ['resid', 'resname', 'atomname', 'distance']
    if oxynotprotein:
        oxy = protein.select_atoms('name O*')
    else:
        oxy = protein.select_atoms('protein and name O*')
    if periodic and (protein.dimensions[:3] > 2).all():
        box = protein.dimensions 
        distances = mda.lib.distances.distance_array(ion.position[np.newaxis, :],
                                                     oxy.positions, box = box)
    else:
        distances = mda.lib.distances.distance_array(ion.position[np.newaxis, :],
                                                     oxy.positions)
    df = pd.DataFrame({'resid': oxy.resids, 'resname': oxy.resnames,
                     'atomname': oxy.names, 'distance': distances[0]},
                     columns=columns)
    df = df[df['distance'] < maxdistance]
    df = df.reset_index()[columns]
    return df

def cume(files, maxdistance=20, binnumber=100, nummols=None):
    """Creates a cumulative histogram of distances of oxygen atoms from an ion.
    :Arguments:
        *files*
            list of locations of files containing distance dataframes
        *maxdistance*
            maximum distance of interest from the ion; default = 20
        *binnumber*
            number of desired bins for cumulative histogram; default = 100
        *nummols*
            number of ions/molecules serving as centers contributing to df; default = None, becomes number of files used

    :Returns:
        *m*
            midpoints of bins
        *cumulative*
            cumulative histogram values
    """
    dataframe = pd.DataFrame()
    x = 0

    for fil in files:
        try:
            f = pd.read_csv(fil, index_col=0)
            dataframe = pd.concat([dataframe, f])
            x += 1
        except:
            with open('failures.out', 'a') as f:
                f.write(fil + '\n')

    if nummols is None:
        nummols = x

    h, e = np.histogram(dataframe[dataframe['distance'] < maxdistance]['distance'], bins=binnumber)
    h = h / float(nummols)
    cumulative = np.cumsum(h)
    m = .5 * (e[:-1] + e[1:])
    return m, cumulative

def gee(dataframes, binnumber=100, nummols=None):
    '''Produces a graph of density as a function of distance
    :Arguments:
        *dataframes*
            list of distance dataframes
        *binnumber*
            number of desired bins for cumulative histogram; default = 100
        *nummols*
            number of ions/molecules serving as centers contributing to df; default = None, becomes number of files used
    :Returns:
        *m*
            midpoints of bins
        *density*
            density histogram values
    '''
    dataframe = pd.concat(dataframes)

    if nummols is None:
        nummols = len(dataframes)

    h, e = np.histogram(dataframe['distance'], bins=binnumber)
    m = .5 * (e[:-1] + e[1:])
    V = 4. / 3 * np.pi * (e[1:] ** 3 - e[:-1] ** 3)

    density = h / V / float(nummols)

    return m, density
