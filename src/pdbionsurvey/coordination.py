# PDB Ion Survey 
# Copyright (c) 2016 Kacey Clark
# Published under the GPL v3
# https://github.com/Becksteinlab/PDB_Ion_Survey/

'''
Functions for analyzing ion coordination in PDB structures
'''

import os.path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import warnings

def en(prot, ion, atomname='O', atomselection='name O* and not name OS', mindistance=.5, maxdistance=20, oxynotprotein=True, periodic=True, pqr=False):
    '''Gives the distances of atoms from an ion.
    :Arguments:
        *prot*
            Treant protein with .pdb or .pqr file
        *ion*
            Atom ion
        *atomname*
            String name of atom; default='O'
        *atomselection*
            String selection for coordinating atoms; default='name O* and not name OS'
        *maxdistance*
            Float maximum distance of interest from the ion; default=20
        *oxynotprotein*
            Boolean true if including oxygens not in the protein; default=True
        *pqr*
            Boolean true if using .pqr file, .pdb else; default=False
    :Returns:
        *df*
            pandas.DataFrame dataframe containing resids, resnames, and atom names for each oxygen in the protein file
    '''
    columns = ['resid', 'resname', 'atomname', 'distance']

    if pqr:
        protein = mda.Universe(prot[prot.name+'.pqr'].abspath)
    else:
        protein = mda.Universe(prot[prot.name+'.pdb'].abspath)

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

    prot['coordination/'+ion.name+'/'+atomname+'/'].make()
    df.to_csv(prot['coordination/'+ion.name+'/'+atomname+'/{}.csv'.format(ion.index)].abspath)
    
    return df

def water_en(prot, atomname='O', atomselection='name O* and not name OS', mindistance=.5, maxdistance=20, oxynotprotein=True, periodic=True, pqr=False):
    '''Gives the distances of atoms from an ion in a waterbox.
    :Arguments:
        *prot*
            Treant protein with .pdb or .pqr file
        *ion*
            Atom ion
        *atomname*
            String name of atom; default='O'
        *atomselection*
            String selection for coordinating atoms; default='name O* and not name OS'
        *maxdistance*
            Float maximum distance of interest from the ion; default=20
        *oxynotprotein*
            Boolean true if including oxygens not in the protein; default=True
        *pqr*
            Boolean true if using .pqr file, .pdb else; default=False
    :Returns:
        *df*
            pandas.DataFrame dataframe containing resids, resnames, and atom names for each oxygen in the protein file
    '''
    columns = ['resid', 'resname', 'atomname', 'distance']

    if pqr:
        protein = mda.Universe(prot[prot.name+'.pqr'].abspath)
    else:
        protein = mda.Universe(prot[prot.name+'.pqr'].abspath)

    if oxynotprotein:
        oxy = protein.select_atoms(atomselection)
    else:
        oxy = protein.select_atoms('protein and '+atomselection)

    origins = protein.select_atoms('resname HOH')

    prot['coordination/WATER/'+atomname+'/'].make()

    dfs = []

    for origin in origins:
        if periodic and (protein.dimensions[:3] > 2).all():
            box = protein.dimensions
            distances = mda.lib.distances.distance_array(origin.position[np.newaxis, :],
                    oxy.positions, box = box)
        else:
            distances = mda.lib.distances.distance_array(origin.position[np.newaxis, :],
                    oxy.positions)

        df = pd.DataFrame({'resid': oxy.resids, 'resname': oxy.resnames,
            'atomname': oxy.names, 'distance': distances[0]}, columns=columns)
        df = df[df['distance'] < maxdistance]
        df = df[df['distance'] > mindistance]
        df = df.reset_index()[columns]

        df.to_csv(prot['coordination/WATER/'+atomname+'/{}.csv'.format(origin.index)].abspath)

    return dfs

def avg_en(bundle, ionname, atomname='O', binsize=.1, nummols=None):
    '''Gives the average distances of atoms from an ion across a bundle.
    :Arguments:
        *bundle*
            Bundle treants
        *ionname*
            String name of ion of interest
        *atomname*
            String name of coordinating atom of interest; defaut='O'
        *binsize*
            Float size of desired bins for histogram; default=.1
        *nummols*
            Float number of ions/molecules serving as centers contributing to df, number of files if None; default=None
    :Returns:
        *ndf*
            pd.DataFrame dataframe containing midpoints of bins and average count over bins
    '''
    frames = []
    for tre in bundle:
        for csv in tre.glob('coordination/'+ionname.upper()+'/'+atomname+'/*.csv'):
            df = pd.read_csv(csv.abspath)
            frames.append(df)

    dataframe = pd.concat(frames)

    if nummols is None:
        nummols = len(frames)

    h, e = np.histogram(dataframe['distance'], bins=np.arange(0, max(dataframe['distance']), binsize))
    m = .5 * (e[:-1] + e[1:])

    n = h / float(nummols)

    ndf = pd.DataFrame({'radius': m, 'count': n}, columns=['radius', 'count'])

    return ndf

def dee(bundle, ionname, atomname='O', binsize=.1, nummols=None):
    '''Gives the number density of atoms around an ion as a function of distance.
    :Arguments:
        *bundle*
            Bundle treants
        *ionname*
            String name of ion of interest
        *atomname*
            String name of coordinating atom of interest; default='O'
        *binsize*
            Float size of desired bins for histogram; default=.1
        *nummols*
            Float number of ions/molecules serving as centers contributing to df, number of files if None; default=None
    :Returns:
        *gdf*
            pd.DataFrame dataframe containing midpoints of bins and number density over bins
    '''
    frames = []
    for tre in bundle:
        for csv in tre.glob('coordination/'+ionname.upper()+'/'+atomname+'/*.csv'):
            df = pd.read_csv(csv.abspath)
            frames.append(df)

    dataframe = pd.concat(frames)

    if nummols is None:
        nummols = len(frames)

    h, e = np.histogram(dataframe['distance'], bins=np.arange(0, max(dataframe['distance']), binsize))
    m = .5 * (e[:-1] + e[1:])
    V = 4. / 3 * np.pi * (e[1:] ** 3 - e[:-1] ** 3)
    density = h / V / float(nummols)

    gdf = pd.DataFrame({'radius': m, 'density': density}, columns=['radius', 'density'])

    return gdf

ATOMNAMES = ['O', 'N', 'S', 'C']

def make_dees(ionname, atomnames=ATOMNAMES, bs=.1, mindistance=True, maxdistance=15, ts=1):
    for atomname in atomnames:
        print('started g '+ionname+' with '+atomname)
        gdf = pdbionsurvey.coordination.gee(b, ionname, atomname=atomname, binsize=bs)
        gdf = gdf[gdf['radius'] < maxdistance]
        print('made g '+ionname+' with '+atomname)
        if not mindistance:
            gdf.to_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
        else:
            mindistance = .5
            gdf['density'] = [gdf['density'][i] if gdf['radius'][i]>mindistance else 0 for i in range(len(gdf['density']))]
            gdf.to_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-withmin.csv')
        print('saved g '+ionname+' with '+atomname)

def closest_oxy_distance(bundle, ion, atom='O', num_oxy=6):
    '''Gives distances of closest oxygens.
    :Arguments:
        *bundle*
            Bundle treants
        *ions*
            String ion name
        *atom*
            String coordinating atom name; default='O'
        *num_oxy*
            Int number of close atoms of interest; default=6
    :Returns:
        *oxys*
            pd.DataFrame dataframe containing distances for the first num_oxy atoms from the ions in ions
    '''
    c = bundle

    z = c[c.tags[ion]]
    oxy = []
    for tre in z:
        for csv in tre.glob('coordination/'+ion.upper()+'/'+atom+'/*.csv'):
            df = pd.read_csv(csv.abspath)
            df = df.sort_values('distance').iloc[:num_oxy]['distance'].values.reshape(1, -1)
            index = '{}_{}'.format(tre.name, csv.name.replace('.csv', ''))
            oxy.append(pd.DataFrame(df, columns=range(1, num_oxy+1), index=[index]))
    oxys = pd.concat(oxy)

    return oxys

def get_peaks(bundle, ionname, mindist=1):
    '''Gives locations of peaks and troughs in g(r)s.
    :Arguments:
        *bundle*
            Bundle treants
        *ionname*
            String name of ion of interest
        *mindist*
            Float minimum distance between peaks and troughs; default=1
    :Returns:
        *m*
            np.array midpoints of bins
        *density*
            np.array density histogram values
        *peaks*
            np.array indices of peak locations
        *mins*
            np.array indices of minimum locations
    '''
    m, density = coordination.gee(bundle, ionname, binnumber=200)
    x = int(round(mindist / (m[1] - m[0])))
    peaks = peakutils.indexes(density, thres=.1, min_dist=x)
    mins = peakutils.indexes(-density, thres=.1, min_dist=x)
    return m, density, peaks, mins

def get_charge(ionname):
    '''Gives locations of peaks and troughs in g(r)s.
    :Arguments:
        *ionname*
            String name of ion of interest
    :Returns:
        *num*
            Int charge of ion
    '''
    if ionname.upper() in ['LI', 'NA', 'K', 'RB', 'CS', 'TL', 'RH', 'AG', 'AU']:
        num = 1
    elif ionname.upper() in ['MG', 'CA', 'SR', 'BA', 'MN', 'CO', 'NI', 'PD', 'PT', 'CU', 'ZN', 'CD', 'HG', 'PB']:
        num = 2
    elif ionname.upper() in ['LA',  'V', 'CR', 'FE', 'RU', 'OS', 'AL', 'GA', 'IN', 'SB']:
        num = 3
    elif ionname.upper() in ['ZR', 'IR']:
        num = 4
    elif ionname.upper() in ['W']:
        num = 6
    elif ionname.upper() in ['F', 'CL', 'BR']:
        num = -1
    elif ionname.upper() in ['IOD', 'I']:
        warnings.warn('Iodide has name I and resname IOD.')
        num = -1
    return num

def set_UkT(prot, ionname, ionselection):
    '''Gives potential energy of protein.
    :Arguments:
        *prot*
            Treant protein
        *ionname*
            String name of ion of interest
        *ionselection*
            String selection of ion
    :Returns:
        *potbyprot*
            List protein name and potential energy
    '''

    ioncharge = get_charge(ionname)
    if not prot[prot.name+'.pqr'].exists:
        return None

    u = mda.Universe(prot[prot.name+'.pdb'].abspath)
    u2 = mda.Universe(prot[prot.name+'.pqr'].abspath)

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

    prot.categories[ionname+str(j)+'_U_kT'] = UkT

    potens.append(poten)
    potbyprot = [prot.name, potens]
    return potbyprot
