"""
Functions for analyzing ion coordination in PDB structures
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import os.path

def en(protein, ion, maxdistance = 20, oxynotprotein = True):
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
    u = protein
    if oxynotprotein:
        oxy = u.select_atoms('name O*')
    else:
        oxy = u.select_atoms('protein and name O*')
    box = u.dimensions
    distances = mda.lib.distances.distance_array(ion.position[np.newaxis, :], 
                                               oxy.positions, box = box)
    oxy_rnames = [atom.resname for atom in oxy]
    oxy_rids = [atom.resid for atom in oxy]
    oxy_names = [atom.name for atom in oxy]    
    df = pd.DataFrame({'resid': oxy_rids, 'resname': oxy_rnames, 
                     'atomname': oxy_names, 'distance': distances[0]},
            columns = ['resid', 'resname', 'atomname', 'distance'])
    df = df[df['distance'] < maxdistance]
    return df

def cume(df, yaxis = 'distance', maxdistance = 20, binnumber = 20, ax = None):
    """Creates a cumulative histogram of distances of oxygen atoms from an ion.
    :Arguments:
        *df*
            `pandas.DataFrame` containing resids, resnames, and atom names
            for each oxygen surrounding the ion
        *yaxis*
            collumn in df to be graphed on the y-axis; default = distance
        *maxdistance*
            maximum distance of interest from the ion; default = 20
        *binnumber*
            number of desired bins for cumulative histogram; default = 20
        *ax*
            axis to plot on; default = None

    :Prints:
        *cumulative histogram*
            number of oxygens within a radius of the ion
        *distances table*
            atom numbers, resids, resnames, atom names, and distances from the ion
            of each oxygen

    :Returns:
        *ax*
            axis used for plotting
    """
    if ax is None:
        fig = plt.figure(figsize = (4,3))
        ax = fig.add_subplot(1,1,1)
    values, base = np.histogram(df[df[yaxis] < maxdistance][yaxis], bins = binnumber)
    cumulative = np.cumsum(values)
    #print df[df['distance'] < maxdistance].sort(columns = 'distance', inplace = False)
    ax.plot(base[:-1], cumulative)
    return ax

def cumin(protein, ions, yaxis = 'distance', maxdistance = 20, oxynotprotein = True, binnumber = 20, ax = None):
    """Creates a cumulative histogram of distances of oxygen atoms from an ion.
    :Arguments:
        *protein*
            protein Universe
        *ions*
            AtomGroup of ions
        *yaxis*
            collumn in df to be graphed on the y-axis; default = distance
        *maxdistance*
            maximum distance of interest from the ion; default = 20
        *oxynotprotein*
            boolean value of whether to include oxygens not in the protein; default = False
        *binnumber*
            number of desired bins for cumulative histogram; default = 20
        *ax*
            axis to plot on; default = None

    :Prints:
        *cumulative histogram*
            number of oxygens within a radius of the ion
        *distances table*
            atom numbers, resids, resnames, atom names, and distances from the ion of each oxygen

    :Returns:
        *ax*
            axis used for plotting
    """
    fig = plt.figure(figsize = (4,3))
    ax = fig.add_subplot(1,1,1)
    for ion in ions:
        cume(en(protein, ion, maxdistance, oxynotprotein), yaxis, maxdistance, binnumber, ax)
    return ax

dataframe = []
proteinids = []

def aggregate(pdbids, path, ionname, maxdistance = 20, oxynotprotein = True):
    """Aggregates dataframes into one dataframe
    :Arguments:
        *pdbids*
            list of PDB codes corresponding to .pdb files to be aggregated (note that .pdb is not to be included in pdbids)
        *path*
            path to the file's directory
        *ionname*
            name of reference ion
        *maxdistance*
            maximum distance of interest from the ion; default = 20
        *oxynotprotein*
            boolean value of whether to include oxygens not in the protein; default = True
    :Returns:
        *dataframe*
            aggregated `pandas.DataFrame` containing resids, resnames, and atom names
            for each oxygen in the protein file
    """
    for x in range(len(pdbids)):
        try:
            u = mda.Universe(os.path.join(path, pdbids[x]) + '.pdb', guess_bonds = False, permissive = False)
            ions = u.select_atoms('not protein and name ' + ionname + '*')
            for i, ion in enumerate(ions):
                dataframe.append(en(u, ion, maxdistance = maxdistance, oxynotprotein = oxynotprotein))
                proteinids.append(pdbids[x] + '_{}'.format(i))
                return dataframe
        except KeyError:
            continue
    return dataframe, proteinids

def aggregategraph(pdbids, path, ionname, yaxis = 'distance', binnumber = 20, maxdistance = 20, oxynotprotein = True):
    """Produces an aggregated graph from multiple dataframes
    :Arguments:
        *pdbids*
            list of PDB codes corresponding to .pdb files to be aggregated (note that .pdb is not to be included in pdbids)
        *path*
            path to the file's directory
        *ionname*
            name of reference ion
        *binnumber*
            number of desired bins for cumulative histogram; default = 20
        *maxdistance*
            maximum distance of interest from the ion; default = 20
        *oxynotprotein*
            boolean value of whether to include oxygens not in the protein; default = True
    :Returns:
        *graph*
            aggregated graph
    """
    dataframe = aggregate(pdbids, path, ionname = ionname, maxdistance = maxdistance, oxynotprotein = oxynotprotein)
    df = pd.concat(dataframe, keys = proteinids, names = ['pdbids'])
    df[yaxis].plot(kind = 'hist', bins = binnumber)
