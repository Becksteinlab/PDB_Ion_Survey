# PDB Ion Survey 
# Copyright (c) 2016 Kacey Clark
# Published under the GPL v3
# https://github.com/Becksteinlab/PDB_Ion_Survey/

"""
Functions for analyzing ion coordination in PDB structures
"""

import os.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import warnings

def en(protein, ion, atomname='O', atomselection='name O* and not name OS', mindistance=.5, maxdistance=20, oxynotprotein=True, periodic=True, sim=None):
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
        *sim*
            sim in which to store df; default=None
    :Returns:
        *df*
        `pandas.DataFrame` containing resids, resnames, and atom names for each oxygen in the protein file
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
    df = df[df['distance'] > mindistance]
    df = df.reset_index()[columns]

    if sim is not None:
        sim['coordination/'+ion.name+'/'+atomname+'/'].make()
        df.to_csv(sim['coordination/'+ion.name+'/'+atomname+'/{}.csv'.format(ion.index)].abspath)
    
    return df

def avg_en(bundle, ionname, atomname='O', binnumber=200, nummols=None):
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
            *n*
                average histogram values
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

    n = h / float(nummols)

    return m, n

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

def closest_oxy_distance(bundle, ion, atom='O', num_oxy=6):
    """Finds distances of closest oxygens.
        :Arguments:
            *bundle*
                bundle of sims
            *ions*
                string ion name
            *atom*
                string coordinating atom name
            *num_oxy*
                number of close oxygens of interest; default = 6
        :Returns:
            *dfs*
                list of DataFrames containing distances for the first num_oxy oxygen atoms from the ions in ions
        """
    c = bundle

    z = c[c.tags[ion]]
    oxy = []
    for sim in z:
        for csv in sim.glob('coordination/'+ion.upper()+'/'+atom+'/*.csv'):
            df = pd.read_csv(csv.abspath)
            df = df.sort_values('distance').iloc[:num_oxy]['distance'].values.reshape(1, -1)
            index = '{}_{}'.format(sim.name, csv.name.replace('.csv', ''))
            oxy.append(pd.DataFrame(df, columns=range(1, num_oxy+1), index=[index]))
    oxys = pd.concat(oxy)

    return oxys

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

def get_charge(ionname):
    if ionname.upper() in ['LI', 'NA', 'K', 'RB', 'CS', 'TL', 'RH', 'AG', 'AU']:
        return 1
    elif ionname.upper() in ['MG', 'CA', 'SR', 'BA', 'MN', 'CO', 'NI', 'PD', 'PT', 'CU', 'ZN', 'CD', 'HG', 'PB']:
        return 2
    elif ionname.upper() in ['LA',  'V', 'CR', 'FE', 'RU', 'OS', 'AL', 'GA', 'IN', 'SB']:
        return 3
    elif ionname.upper() in ['ZR', 'IR']:
        return 4
    elif ionname.upper() in ['W']:
        return 6
    elif ionname.upper() in ['F', 'CL', 'BR']:
        return -1
    elif ionname.upper() in ['IOD', 'I']:
        warnings.warn('Iodide has name I and resname IOD.')
        return -1

def set_UkT(sim, ionname, ionselection):
    '''
    :Arguments:
        *sim*
            sim
        *ionname*
            string name of ion
        *ionselection*
            string selection of ion
    '''

    ioncharge = get_charge(ionname)
    if not sim[sim.name+'.pqr'].exists:
        return None

    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    u2 = mda.Universe(sim[sim.name+'.pqr'].abspath)

    ions = u.select_atoms(ionselection)
    atoms = u2.atoms

    potens = []

    if len(ions) == 0:
        return None

    for j, ion in enumerate(ions):
        columns = ['resid', 'resname', 'atomname', 'charge', 'distance', 'coulomb potential']
        box = u.dimensions
        distances = mda.lib.distances.distance_array(ion.position[np.newaxis, :],
                    atoms.positions, box = box)
        potentials = []
        for i, atom in enumerate(atoms):
            potentials.append(atom.charge * ioncharge * (1.60217662e-19) * (1.60217662e-19) * 8.9875517873681764e9 / (distances[0,i] * 10**-10))
        poten = pd.DataFrame({'resid': atoms.resids, 'resname': atoms.resnames, 'atomname': atoms.names, 'charge': atoms.charges, 'distance': distances[0], 'coulomb potential': potentials}, columns=columns)

    potenergy = sum(poten['coulomb potential'])

    kb = 1.38064852e-23

    kT = kb * 300

    UkT = potenergy/kT

    sim.categories[ionname+str(j)+'_U_kT'] = UkT

    potens.append(poten)
    return [sim.name, potens]
