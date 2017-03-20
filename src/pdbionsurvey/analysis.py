# PDB Ion Survey 
# Copyright (c) 2016 Kacey Clark
# Published under the GPL v3
# https://github.com/Becksteinlab/PDB_Ion_Survey/

"""
Functions for creating sims and analyzing coordination data.
"""

from __future__ import absolute_import

import os

import mdsynthesis as mds
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import peakutils

from . import coordination

def make_sims(path, pdbfiles):
    """Makes a Tree of sims in a path.
    :Arguments:
        *path*
            path to sims
        *pdbfiles*
            list of paths to .pdb files

    :Returns:
        *sims*
            Tree of sims
        *pdbfiles*
            View of .pdb files
    """
    sims = mds.Bundle()
    for fil in pdbfiles:
        mds.Sim(path.abspath+fil.name[:4])
        sims.add(path.abspath+fil.name[:4])
        with open(sims[fil.name[:4]][0].abspat, 'w') as f:
            f.write(fil.read())
            f.close()
    return sims

def define_universe(sims, pdbfiles):
    """Makes universes in sims.
    :Arguments:
        *sims*
            Tree of sims
        *pdbfiles*
            list of paths to .pdb files
    """
    for fil in pdbfiles:
        try:
            sims[fil.name[:4]][0].universe = mda.Universe(fil.abspath)
        except:
            with open('failures.out', 'a') as f:
                f.write(fil.name + '\n')

def sim_labeling(bundle, ionname=None, project_tags=['pdbionsurvey', 'pdbsurvey']):
    """Adds tags and categories to sims.
    :Arguments:
        *bundle*
            Bundle object
        *ionname*
            name of ion for tagging; default = None
        *project_tags*
            list of tags to attach to sims; default = ['pdbionsurvey', 'pdbsurvey']
    """
    for tag in project_tags:
        bundle.tags.add(tag)

    if ionname: 
        bundle.tags.add(ionname)
        bundle.tags.add(ionname.upper())
        bundle.tags.add(ionname.lower())

    for sim in bundle:
        with open(sim.glob('*.pdb')[0].abspath, "r") as f:
           searchlines = f.readlines()
        found = False

        for i, line in enumerate(searchlines):
            if ("REMARK" in line) and ("RESOLUTION." in line):
                try:
                    things = line.split()
                    resolutionindex = things.index('RESOLUTION.')
                    value = things[resolutionindex + 1]
                    value = float(value)
                    found = True
                    break
                except (IndexError, KeyError):
                    continue
                except ValueError:
                    continue

        if found:
            sim.categories['resolution'] = value
        else:
            sim.categories['resolution'] = 'N/A'

        if sim.universe.select_atoms('resname HOH'):
            value = True
        else:
            value = False

        sim.categories['has_water'] = value

def closest_oxy_distance(bundle, ions, resolutions, cume = True, num_oxy = 6):
    """Finds distances of closest oxygens.
    :Arguments:
        *bundle*
            bundle of sims
        *ions*
            list of ion names
        *resolutions*
            list of resolution numbers
        *cume*
            boolean value of whether to sort by resolution cumulatively; default = True
        *num_oxy*
            number of close oxygens of interest; default = 6
    :Returns:
        *dfs*
            list of DataFrames containing distances for the first num_oxy oxygen atoms from the ions in ions
    """
    c = bundle
    dfs = []
    for res in resolutions:
        if not cume:
            c = c[[r > (res - .5) for r in c.categories['resolution']]]
        c = c[[r <= res for r in c.categories['resolution']]]

        for ion in ions:
            z = c[c.tags[ion]]
            oxy = []
            for sim in z:
                for csv in sim.glob('coordination/'+ion.upper()+'/*.csv'):
                    df = pd.read_csv(csv.abspath)
                    df = df.sort_values('distance').iloc[:num_oxy]['distance'].values.reshape(1, -1)
                    index = '{}_{}'.format(sim.name, csv.name.replace('.csv', ''))
                    oxy.append(pd.DataFrame(df, columns=range(1, num_oxy+1), index=[index]))
            oxys = pd.concat(oxy)
            dfs.append(oxys)
    return dfs

def graph_closest_oxy_distances(dfs, ax=None, cume=True, axlim=(1, 6), binsize = .2, save=False):
    """Creates a neat plot of closest oxygen distance data.
    :Arguments:
        *dfs*
            pandas.DataFrame` containing distances for the first num_oxy oxygen atoms from the ions in ions
        *ax*
            axes object; default = None
        *cume*
            boolean value of whether to sort by resolution cumulatively; default = True
        *axlim*
            minimum and maximimum distances from ion of interest; default = (1, 6)
         *binsize*
            bin width of distances from ion; default = .2
        *save*
            boolean value of whether to save graph; default = False
    :Returns:
        *ax*
            axes object
    """
    if not ax:
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(1,1,1)
    ax.set_xlim(axlim)

    bins = np.arange(0, 8, binsize)

    for i in range(len(dfs.columns)):
        h, e = np.histogram(dfs[i+1], bins=bins)
        m = .5 * (e[:-1] + e[1:])
        ax.plot(m, h, label='oxy #' + str(i+1), lw=2)

    ax.set_xlabel('Distance ($\AA$)')
    ax.set_ylabel('Frequency')
    ax.legend(fontsize='medium')

    if save and cume:
        ax.figure.savefig(ion + "+ Distance of Oxygens Histogram (res from " + (res - .5) + " to " + res + ").pdf")
    elif save:
        ax.figure.savefig(ion + "+ Distance of Oxygens Histogram (res to " + res + ").pdf")
    return ax

def get_peaks(bundle, ionname, mindist=1):
    """Calculates location of peaks and troughs in g(r)s.
    :Arguments:
        *bundle*
            bundle of sims
        *ionname*
            name of ion of interest
        *mindist*
            minimum distance between peaks and troughs; default = 1
    :Returns:
        *m*
            midpoints of bins
        *density*
            density histogram values
        *peaks*
            indexes of peak locations
        *mins*
            indexes of minimum locations
    """
    m, density = coordination.gee(bundle, ionname, binnumber=200)
    x = int(round(mindist / (m[1] - m[0])))
    peaks = peakutils.indexes(density, thres=.1, min_dist=x)
    mins = peakutils.indexes(-density, thres=.1, min_dist=x)
    return m, density, peaks, mins
