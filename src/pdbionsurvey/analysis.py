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

IONNAMES = ['NA', 'LI', 'K', 'TL', 'CL', 'ZN', 'CA']

def make_sims(path, pdbid):
    """Makes a Tree of sims in a path.
    :Arguments:
        *path*
            path to sims
        *pdbid*
            PDB id code
    """
    sim = mds.Sim(path.abspath+pdbid)
    i = mmtf.fetch(pdbid)
    u = mda.Universe(i)
    u.atoms.write(sim.abspath+pdbid+'.pdb', bonds=None)
    sim.universe = u
    sim.categories['resolution'] = i.resolution

def sim_labeling(sim, ionnames=IONNAMES, project_tags=['pdbionsurvey', 'pdbsurvey'], ligands=True):
    """Adds tags and categories to sims.
    :Arguments:
        *sim*
            sim object containing a .pdb file
        *ionnames*
            names of ions for tagging; default = ['Na', 'Li', 'K', 'Tl', 'Cl']
        *project_tags*
            list of tags to attach to sims; default = ['pdbionsurvey', 'pdbsurvey']
        *ligands*
            boolean value of whether to look for nonamino acid resnames (includes ions); default = True
    """
    for tag in project_tags:
        sim.tags.add(tag)

    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    for ionname in ionnames:
        ions = u.select_atoms('name '+ionname.upper()+ ' and resname '+ionname.upper())
        if len(ions) != 0:
            sim.tags.add(ionname)
            sim.tags.add(ionname.upper())
            sim.tags.add(ionname.lower())
            sim.categories['num_'+ionname.upper()] = len(ions)

    if ligands:
        if sim[sim.name+'.pqr'].exists:
            try:
                v = mda.Universe(sim[sim.name+'.pqr'].abspath)
                ligs = set(u.atoms.resnames) - set(v.atoms.resnames)
                sim.tags.remove('ligand')
                ligs = list(ligs - ligs.intersection(IONNAMES))
                if len(ligs) != 0:
                    sim.tags.add('ligand')
                for lig in ligs:
                    sim.tags.add(lig)
                sim.tags.remove('Failed to construct topology from .pqr file')
            except ValueError:
                sim.tags.add('Failed to construct topology from .pqr file')

    if not sim[sim.name+'.pqr'].exists:
        sim.tags.add('no_pqr')

    if u.select_atoms('resname HOH'):
        value = True
    else:
        value = False

    sim.categories['has_water'] = value

    sim.categories['smallest_dimension'] = float(min(u.dimensions[:3]))
    if not (u.dimensions[:3] > 2).all():
        sim.tags.add('funky_dimensions')

def pdb2pqr(sim, pdb2pqrpath='/nfs/packages/opt/Linux_x86_64/pdb2pqr/2.1.1/pdb2pqr.py'):
    '''
    :Arguments:
        *sim*
            sim
        *pdb2pqrpath*
            string path to pdb2pqr
    '''
    if not sim[sim.name+'.pqr'].exists:
        try:
            return os.system(pdb2pqrpath+' --ff=charmm --whitespace {} {}'.format(sim.relpath+sim.name+'.pdb', sim.relpath+sim.name+'.pqr'))
        except:
            sim.tags.add('no_pqr')
    else:
        sim.tags.remove('no_pqr')

def pdb2mol2(sim, lig, babelpath='/usr/bin/babel'):
    '''
    :Arguments:
        *sim*
            sim
        *lig*
            string ligand name
        *babelpath*
            string path to babel
    '''
    return os.system(babelpath+' {} {}'.format(sim['ligands/'+lig.upper()+'.pdb'], sim['ligands/'+lig.upper()+'.mol2']))

def ligsolution(sim):
    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    v = mda.Universe(sim[sim.name+'.pqr'].abspath)
    ligs = set(u.atoms.resnames) - set(v.atoms.resnames)
    ligs = list(ligs - ligs.intersection(IONNAMES))
    for lig in ligs:
        ligatoms = u.select_atoms('resname '+lig)
        sim['ligands/'].make()
        ligatoms.write(sim['ligands/'+lig.upper()+'.pdb'].abspath)
        pdb2mol2(sim, lig)

def closest_oxy_distance(bundle, ion, num_oxy=6):
    """Finds distances of closest oxygens.
    :Arguments:
        *bundle*
            bundle of sims
        *ions*
            ion name
        *num_oxy*
            number of close oxygens of interest; default = 6
    :Returns:
        *df*
            DataFrame containing distances for the first num_oxy oxygen atoms
    """
    c = bundle
    z = c[c.tags[ion]]
    oxy = []
    for sim in z:
        for csv in sim.glob('coordination/'+ion.upper()+'/O/*.csv'):
            df = pd.read_csv(csv.abspath)
            df = df.sort_values('distance').iloc[:num_oxy]['distance'].values.reshape(1, -1)
            index = '{}_{}'.format(sim.name, csv.name.replace('.csv', '')
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

def set_UkT(sim, ionname, ionselection, ioncharge=1):
    '''
    :Arguments:
        *sim*
            sim
        *ionname*
            string name of ion
        *ionselection*
            string selection of ion
        *ioncharge*
            int or float charge of ion in e-
    '''

    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    u2 = mda.Universe(sim[sim.name+'.pqr'].abspath)

    ions = u.select_atoms(ionselection)
    atoms = u2.atoms

    potens = []
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
