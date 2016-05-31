"""
Functions for analyzing ion coordination in PDB structures
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import os.path

def en(protein, ion, maxdistance=20, oxynotprotein=True):
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
            columns=columns)
    df = df[df['distance'] < maxdistance]
    df = df.reset_index()[columns]
    return df

def cume(df, yaxis='distance', maxdistance=20, binnumber=20, nummols=1, ax=None):
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
            number of desired bins for cumulative histogram; default = 20        *ax*
            axis to plot on; default = None
        *nummols*
            number of proteins/molecules contributing to df

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

    values = values/float(nummols)
    cumulative = np.cumsum(values)
    #print df[df['distance'] < maxdistance].sort(columns = 'distance', inplace = False)
    ax.plot(base[:-1], cumulative)
    return ax

def cumin(protein, ions, yaxis='distance', maxdistance=20, oxynotprotein=True,
          binnumber=20, nummols=1, ax=None):
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
        cume(en(protein=protein, ion=ion, maxdistance=maxdistance, oxynotprotein=oxynotprotein), yaxis=yaxis, maxdistance=maxdistance, binnumber=binnumber, nummols=nummols, ax=ax)
    return ax

def gee(files, filename):
    '''Produces a graph of density as a function of distance
    :Arguments:
        *files*
            list of dataframe file locations
        *filename*
            desired name of output file
    :Returns:
        *graph*
            graph of g(r)
    '''
    dataframe2 = []
    for x in range(len(files)):
        try:
            f = pd.read_csv(pat + files[x], index_col=0)
            pdbid = f[0:4]
            dataframe2.append(f)
        except:
            with open('failures.out', 'a') as f:
                f.write(pdbfile.name + '\n')

    df = pd.concat(dataframe2, keys = pdbid, names = ['pdbids'])

    h, e = np.histogram(df['distance'], bins = 100)
    m = .5 * (e[:1] + e[1:])
    V = 4 / 3 * np.pi * (e[1:] ** 3 - e[:1] ** 3)

    density = h / V

    ax = plt.subplot(1,1,1)
    ax.plot(m, density)

def aggregate(pdbids, path, ionname, maxdistance=20, oxynotprotein=True):
    """Aggregates dataframes into one dataframe
    :Arguments:
        *pdbids*
            list of PDB codes corresponding to .pdb files to be aggregated
            (note that .pdb is not to be included in pdbids) 
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
    dataframe = []

    for x in range(len(pdbids)):
        try:
            u = mda.Universe(os.path.join(path, pdbids[x]) + '.pdb', guess_bonds = False, permissive = False)
            ions = u.select_atoms('not protein and name ' + ionname + '*')
            for i, ion in enumerate(ions):
                dataframe.append(en(u, ion, maxdistance = maxdistance, oxynotprotein = oxynotprotein))
                proteinids.append(pdbids[x] + '_{}'.format(i))
                return dataframe
        except IOError:
            continue
    return dataframe, proteinids

def aggregategraph(pdbids, path, ionname, yaxis='distance',
                   binnumber=20, maxdistance=20, oxynotprotein=True):
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
    dataframe = aggregate(pdbids, path, ionname=ionname,
                          maxdistance=maxdistance, oxynotprotein=oxynotprotein)
    df = pd.concat(dataframe, keys=proteinids, names=['pdbids'])
    df[yaxis].plot(kind='hist', bins=binnumber)
