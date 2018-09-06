import os

import mdsynthesis as mds

import pdbionsurvey.analysis as analysis

analysis.make_sims??

os.mkdir('sims/')
path = 'sims/'
analysis.make_sims('2JLN', path)
sim = mds.Sim(path+'2JLN')
analysis.pdb2pqr(sim)
analysis.ligsolution(sim)
analysis.allligs(sim)
analysis.pdb2pqrcomplete(sim)
