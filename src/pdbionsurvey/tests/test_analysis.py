import os
import pytest

import mdsynthesis as mds

import pdbionsurvey.analysis as analysis

@pytest.fixture
def sim(tmpdir):
    p = tmpdir.mkdir('sims/')
    path = str(p)
    analysis.make_sims('1AK0', path=path)
    sim = mds.Sim(os.path.join(path, '1AK0'))
    return sim

def test_simmaking(sim):
    assert(sim[sim.name+'.pdb'].exists)

def test_pqrmaking(sim):
    analysis.pdb2pqr(sim)
    assert(sim[sim.name+'.pqr'].exists)

def test_ligmaking(sim):
    analysis.pdb2pqr(sim)
    analysis.ligsolution(sim)
    assert(sim['ligands/'].exists)

def test_allligs(sim):
    analysis.pdb2pqr(sim)
    analysis.ligsolution(sim)
    analysis.allligs(sim)
    assert(sim['ligands/all_ligands.pdb'].exists)
