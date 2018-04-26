# PDB Ion Survey 
# Copyright (c) 2016 Kacey Clark
# Published under the GPL v3
# https://github.com/Becksteinlab/PDB_Ion_Survey/

"""
Functions for analyzing ion coordination in PDB structures
"""

from __future__ import absolute_import

import os.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

def en(protein, ion, atomname='O', atomselection='name O* and not name OS', maxdistance=20, oxynotprotein=True, periodic=True):
    """Gives the distances of oxygen atoms from an ion.
    :Arguments:
    *protein*
        protein Universe
    *ion*
        ion Atom
    *atomname*
        string name of atom
    *atomselection*
        string selection for coordinating atoms
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
        oxy = protein.select_atoms(atomselection)
    else:
        oxy = protein.select_atoms('protein and '+atomselection)
    if periodic and (protein.dimensions[:3] > 2).all():
        box = protein.dimensions 
        distances = mda.lib.distances.distance_array(ion.position[np.newaxis, :],
                    oxy.positions, box = box)
    else:
        distances = mda.lib.distances.distance_array(ion.position[np.newaxis, :],
                    oxy.positions)
    df = pd.DataFrame({'resid': oxy.resids, 'resname': oxy.resnames,
            'atomname': oxy.names, 'distance': distances[0]}, columns=columns)
    df = df[df['distance'] < maxdistance]
    df = df.reset_index()[columns]

    df.to_csv(sim['coordination/'+ion+'/'+atomname+'/{}.csv'.format(ion.index)].abspath)

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

def gee(bundle, ionname, atomname='O', binnumber=200, nummols=None):
    '''Produces a graph of density as a function of distance
    :Arguments:
    *bundle*
    bundle of sims
    *ionname*
    name of ion of interest
    *atomname*
    name of coordinating atom of interest
    *binnumber*
    number of desired bins for cumulative histogram; default = 200
    *nummols*
    number of ions/molecules serving as centers contributing to df; default = None, becomes number of files used
    :Returns:
    *m*
    midpoints of bins
    *density*
    density histogram values
    '''
    frames = []
    for sim in bundle:
        for csv in sim.glob('coordination/'+ionname.upper()+'/'+atomname+'/*.csv'):
            df = pd.read_csv(csv.abspath)
            frames.append(df)

    dataframe = pd.concat(frames)

    if nummols is None:
        nummols = len(frames)

    h, e = np.histogram(dataframe['distance'], bins=binnumber)
    m = .5 * (e[:-1] + e[1:])
    V = 4. / 3 * np.pi * (e[1:] ** 3 - e[:-1] ** 3)

    density = h / V / float(nummols)

    return m, density
