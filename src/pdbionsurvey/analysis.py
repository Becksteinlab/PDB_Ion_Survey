# PDB Ion Survey 
# Copyright (c) 2016 Kacey Clark
# Published under the GPL v3
# https://github.com/Becksteinlab/PDB_Ion_Survey/

'''
Functions for creating sims and analyzing coordination data.
'''

from __future__ import absolute_import

import os
import subprocess

import datreant as dtr
import MDAnalysis as mda
import numpy as np
import pandas as pd
import peakutils
import mmtf

from . import coordination

IONNAMES = ['NA', 'MG', 'K', 'CA', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'PD', 'AG', 'CD', 'IR', 'PT', 'AU', 'HG', 'LA', 'PB', 'TL', 'LI', 'BA', 'RB', 'CS', 'SR', 'CL', 'IOD', 'F', 'BR']

def make_treants(pdbid, path='proteins/'):
    '''Makes a Tree of treants in a path.
    :Arguments:
        *pdbid*
            String pdb id code
        *path*
            String path to treants; default='proteins/'
    '''
    tree = dtr.Tree(protdir.abspath+pdbid)
    prot = mmtf.fetch(pdbid)
    u = mda.Universe(prot)
    for i, mod in enumerate(u.models):
        tre = dtr.Treant(tree.abspath+pdbid+str(i))
        mod.atoms.write(tre.abspath+pdbid+str(i)+'.pdb', bonds=None)
        if prot.resolution:
            tre.categories['resolution'] = prot.resolution
        else:
            tre.categories['resolution'] = 'N/A'
        tre.categories['total_models'] = len(u.models)
        tre.tags.add(pdbid)
        tre.categories['pdbid'] = pdbid
        tre.categories['model_num'] = i

def sim_labeling(sim, ionnames=IONNAMES, project_tags=['pdbionsurvey', 'pdbsurvey'], ligands=True):
    '''Adds tags and categories to sims.
    :Arguments:
        *sim*
            Sim sim containing  a .pdb file
        *ionnames*
            List names of ions for tagging; default=IONNAMES
        *project_tags*
            List tags to attach to sims; default=['pdbionsurvey', 'pdbsurvey']
        *ligands*
            Boolean true if looking for nonamino acid resnames (includes ions); default=True
    :Returns:
    '''
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

def pdb2pqr(sim, pdb2pqrloc='/nfs/packages/opt/Linux_x86_64/pdb2pqr/2.1.1/pdb2pqr.py'):
    '''Converts .pdb file to .pqr file.
    :Arguments:
        *sim*
            Sim sim containing .pdb file
        *pdb2pqrloc*
            String location of pdb2pqr; default='/nfs/packages/opt/Linux_x86_64/pdb2pqr/2.1.1/pdb2pqr.py'
    :Returns:
    '''
    return subprocess.run(pdb2pqrloc +' --ff=charmm --whitespace {} {}'.format(sim.relpath+sim.name+'.pdb', sim.relpath+sim.name+'.pqr'))

def pdb2mol2(sim, lig, babelloc='/usr/bin/babel'):
    '''Converts .pdb file to .mol2 file.
    :Arguments:
        *sim*
            Sim sim containing .pdb file
        *babelloc*
            String location of babel; default='/usr/bin/babel'
    :Returns:
    '''

    return subprocess.run(babelloc+' {} {}'.format(sim['ligands/'+lig.upper()+'.pdb'], sim['ligands/'+lig.upper()+'.mol2']))

def ligsolution(sim):
    '''Converts ligs in .pdb file to .mol2 files.
    :Arguments:
        *sim*
            Sim sim containing .pdb and .pqr files
    :Returns:
    '''

    if not sim[sim.name+'.pqr'].exists:
        return None
    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    v = mda.Universe(sim[sim.name+'.pqr'].abspath)
    ligs = set(u.atoms.resnames) - set(v.atoms.resnames)
    ligs = list(ligs - ligs.intersection(set(IONNAMES)))
    for lig in ligs:
        ligatoms = u.select_atoms('resname '+lig)
        sim['ligands/'].make()
        sim.tags.add(lig)
        ligatoms.write(sim['ligands/'+lig.upper()+'.pdb'].abspath)
        pdb2mol2(sim, lig)

def allligs(sim, allname='all_ligands.pdb'):
    '''Puts all ligands in .pdb file.
    :Arguments:
        *sim*
            Sim sim containing ligand files
        *pdb2pqrloc*
            String desired filename
    :Returns:
        '''
    u = mda.Universe(sim[sim.name+'.pdb'].abspath)
    v = None
    if sim['ligands/'].exists:
        for leaf in sim['ligands/'].leaves():
            if leaf.name == allname:
                continue
            lig = u.select_atoms('resname '+leaf.name[:-4])
            if v:
                v = v + (lig - v)
            else:
                v = lig
        v.write(sim[os.path.join('ligands', allname)].abspath)

def pdb2pqrcomplete(sim, pdb2pqrloc='/nfs/packages/opt/Linux_x86_64/pdb2pqr/2.1.1/pdb2pqr.py'):
    '''Converts .pdb file to .pqr file with ligands.
    :Arguments:
        *sim*
            Sim sim containing .pdb file
        *pdb2pqrloc*
            String location of pdb2pqr; default='/nfs/packages/opt/Linux_x86_64/pdb2pqr/2.1.1/pdb2pqr.py'
    :Returns:
    '''
    return subprocess.run(pdb2pqrloc+' --ff=charmm --ligand {} --whitespace {} {}'.format(sim['ligands/all_ligands.mol2'].relpath, sim[sim.name+'-complete.pdb'].relpath, sim[sim.name+'-withligs.pqr'].relpath))
