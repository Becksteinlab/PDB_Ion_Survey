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

IONNAMES = ['NA', 'MG', 'K', 'CA', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'PD', 'AG', 'CD', 'IR', 'PT', 'AU', 'HG', 'LA', 'PB', 'TL', 'LI', 'BA', 'RB', 'CS', 'SR', 'CL', 'IOD', 'F', 'BR']

def make_sims(path='sims/', pdbid):
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
    try:
        sim.categories['resolution'] = i.resolution
    except:
        sim.categories['resolution'] = 'N/A'

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

    if not (u.dimensions[:3] > 2).all():
        sim.tags.add('funky_dimensions')

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

def pdb2pqr(sim):
    return os.system('/nfs/packages/opt/Linux_x86_64/pdb2pqr/2.1.1/pdb2pqr.py --ff=charmm --whitespace {} {}'.format(sim.relpath+sim.name+'.pdb', sim.relpath+sim.name+'.pqr'))

def pdb2mol2(lig):
    return os.system('/usr/bin/babel {} {}'.format(sim['ligands/'+lig.upper()+'.pdb'], sim['ligands/'+lig.upper()+'.mol2']))

def ligsolution(sim):
    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    v = mda.Universe(sim[sim.name+'.pqr'].abspath)
    ligs = set(u.atoms.resnames) - set(v.atoms.resnames)
    ligs = list(ligs - ligs.intersection(set(IONNAMES)))
    for lig in ligs:
        ligatoms = u.select_atoms('resname '+lig)
        sim['ligands/'].make()
        sim.tags.add(lig)
        ligatoms.write(sim['ligands/'+lig.upper()+'.pdb'].abspath)
        pdb2mol2(lig)

def allligs(sim):
    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    v = None
    if sim['ligands/'].exists:
        for leaf in sim['ligands/'].leaves():
            lig = u.select_atoms('resname '+leaf.name[:-4])
            if v:
                v = v + (lig - v)
            else:
                v = lig
        v.write(sim['ligands/all_ligands.pdb'].abspath)

def complete(filename, outname):
    env = environ()
    # topology parameters for building from internal coordinates
    # modified Charmm22
    env.libs.topology.read('$(LIB)/top_heav.lib')
    env.libs.parameters.read('$(LIB)/par.lib')

    mdl = modeller.scripts.complete_pdb(env, filename)
    mdl.write(outname)
    return outname

def completesim(sim):
    complete(sim[sim.name+'.pdb'].abspath, sim[sim.name+'-complete.pdb'].abspath)

def pdb2pqrcomplete(sim):
    return os.system('/nfs/packages/opt/Linux_x86_64/pdb2pqr/2.1.1/pdb2pqr.py --ff=charmm --ligand {} --whitespace {} {}'.format(sim['ligands/all_ligands.mol2'].relpath, sim[sim.name+'-complete.pdb'].relpath, sim[sim.name+'-withligs.pqr'].relpath))
