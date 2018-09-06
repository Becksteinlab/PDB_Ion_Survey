import os
import pytest

import mdsynthesis as mds

import pdbionsurvey.analysis as analysis

def test_ligmaking(tmpdir):
    p = tmpdir.mkdir('sims/')
    path = str(p)
    analysis.make_sims('1AK0', path)
    sim = mds.Sim(path+'1AK0')
    analysis.pdb2pqr(sim)
    assert(sim[sim.name+'.pqr'].exists)
    analysis.ligsolution(sim)
    assert(sim['ligands/'].exists)
    analysis.allligs(sim)
    assert(sim['ligands/all_ligands.pdb'].exists)
