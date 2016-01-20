"""
Functions for analyzing ion coordination in PDB structures
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

def gee(protein, ion, maxdistance=20, oxynotprotein=False):
    """Gives the distances of oxygen atoms from an ion.

    :Arguments:
        *protein*
            protein Universe
        *ion*
            ion Atom
        *maxdistance*
            maximum distance of interest from the ion; default=20
        *oxynotprotein*
            boolean value of whether to include oxygens not in the protein; default=False

    :Returns:
        *df*
            `pandas.DataFrame` containing resids, resnames, and atom names
            for each oxygen in the protein file
    """
    u=protein
    if oxynotprotein:
        oxy=u.select_atoms('name O*')
    else:
        oxy=u.select_atoms('protein and name O*')
    d=oxy.positions-ion.position
    distance=(np.sum(d*d, axis=1))**.5
    distances=list(distance)
    oxy_rnames=[atom.resname for atom in oxy]
    oxy_rids=[atom.resid for atom in oxy]
    oxy_names=[atom.name for atom in oxy]
    df=pd.DataFrame({'resid': oxy_rids, 'resname': oxy_rnames, 'atomname': oxy_names, 'distance': distances},
            columns=['resid', 'resname', 'atomname', 'distance'])
    df=df[df['distance'] < maxdistance]
    return df

def ofr(df, maxdistance=20, binnumber=20, ax=None):
    """Creates a cumulative histogram of distances of oxygen atoms from an ion.

    :Arguments:
        *df*
            `pandas.DataFrame` containing resids, resnames, and atom names
            for each oxygen surrounding the ion
        *maxdistance*
            maximum distance of interest from the ion; default=20
        *binnumber*
            number of desired bins for cumulative histogram; default=20
        *ax*
            axis to plot on; default=None

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
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(1,1,1)
    values, base=np.histogram(df[df['distance'] < maxdistance]['distance'], bins=binnumber)
    cumulative=np.cumsum(values)
    #print df[df['distance'] < maxdistance].sort(columns='distance', inplace=False)
    ax.plot(base[:-1], cumulative)
    return ax

def gofr(protein, ions, maxdistance=20, oxynotprotein=False, binnumber=20, ax=None):
    """Creates a cumulative histogram of distances of oxygen atoms from an ion.

    :Arguments:
        *protein*
            protein Universe
        *ions*
            AtomGroup of ions
        *maxdistance*
            maximum distance of interest from the ion; default=20
       *oxynotprotein*
            boolean value of whether to include oxygens not in the protein; default=False
        *binnumber*
            number of desired bins for cumulative histogram; default=20
        *ax*
            axis to plot on; default=None

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
    fig=plt.figure(figsize=(4,3))
    ax=fig.add_subplot(1,1,1)
    for ion in ions:
        ofr(gee(protein, ion, maxdistance, oxynotprotein), maxdistance, binnumber, ax)
    return ax
