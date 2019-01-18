import MDAnalysis as mda
import mdsynthesis as mds
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import seaborn as sb
from os import path
import os
import shutil
import peakutils
from glob import glob

#table 1 in Peters paper
newcations = ['NA', 'MG', 'K', 'CA', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'PD', 'AG', 'CD', 'IR', 'PT', 'AU', 'HG']
othercations = ['LA', 'PB', 'TL', 'LI', 'BA', 'RB', 'CS', 'SR']
anions = ['CL', 'I', 'F', 'BR']

IONNAMES = list(set(newcations).union(set(othercations), set(anions)))
atoms = ['O', 'N', 'S']

def kolmogorov_smirnov(g1, g2, ylabel='count', cumulative=False):
    assert(np.equal(np.array(g1['radius']), np.array(g2['radius'])).all())
    if not cumulative:
        g1d = np.cumsum(np.array(g1[ylabel]))
        g2d = np.cumsum(np.array(g2[ylabel]))
    else:
        g1d = np.array(g1[ylabel])
        g2d = np.array(g2[ylabel])
    g1d = g1['count']/max(g1['count'])
    g2d = g2['count']/max(g2['count'])
    diff = abs(g2d - g1d)
    return max(diff)

def cramer_vonmises(g1, g2, ylabel='count', cumulative=False):
    assert(np.equal(np.array(g1['radius']), np.array(g2['radius'])).all())
    if not cumulative:
        g1d = np.cumsum(np.array(g1[ylabel]))
        g2d = np.cumsum(np.array(g2[ylabel]))
    else:
        g1d = np.array(g1[ylabel])
        g2d = np.array(g2[ylabel])
    diff = abs(g2d - g1d)
    squarediff = np.square(diff)
    parts = squarediff*(g1d-np.insert(g1d[:-1], 0, [0]))
    wsquared = np.sum(parts)
    return wsquared

def c(alpha):
    c = np.sqrt(-1/2*np.log(alpha))
    return c

def kolmogorov_smirnov_test(g1, g2, ylabel='count', cumulative=False, alpha=.001):
    ks = kolmogorov_smirnov(g1, g2, ylabel='count', cumulative=False)
    n1 = len(g1)
    n2 = len(g2)
    calpha = c(alpha)
    std = calpha * np.sqrt((n1+n2)/(n1*n2))
    if ks > std:
        print('Distinct with ', alpha, 'confidence.')
        return True
    else:
        print('Not distinct with', alpha, 'confidence.')
        return False

def ks_test_comparison(ionnames=IONNAMES, alpha=1):
    ksstats = np.zeros((len(IONNAMES)*3, len(IONNAMES)*3))
    for k in range(len(IONNAMES)):
        for l in range(len(IONNAMES)):
            for i in range(3):
                for j in range(3):
                    atoms = ['O', 'N', 'S']
                    ionname1 = IONNAMES[k]
                    ionname2 = IONNAMES[l]
                    atomname1 = atoms[i]
                    atomname2 = atoms[j]
                    g1 = pd.read_csv('gofrs/n-'+ionname1+'-'+atomname1+'.csv')
                    g1 = pd.DataFrame({'radius': g1['radius'], 'count': np.cumsum(g1['count'])/max(np.cumsum(g1['count']))})
                    g2 = pd.read_csv('gofrs/n-'+ionname2+'-'+atomname2+'.csv')
                    
                    g2 = pd.DataFrame({'radius': g2['radius'], 'count': np.cumsum(g2['count'])/max(np.cumsum(g2['count']))})

                    kst = kolmogorov_smirnov_test(g1, g2, ylabel='count', alpha=alpha)
                    ksstats[k*3+i][l*3+j] = kst

    u = []
    for i in IONNAMES:
        for j in ['O', 'N', 'S']:
            u.append(i+'_'+j)

    dictionary = dict(zip(u, ksstats))
    dictionary['ion/atom'] = u
    df = pd.DataFrame(dictionary, columns=['ion/atom']+u)
    return df

def ks_test_multi(alphas=[.5, .1, .01, .001])
    dfs = []
    for a in alphas:
        df = ks_test_comparison(alpha=a)
        dfs.append(df)
    bigdf = pd.concat(dfs,keys=alphas,axis=0)
    return bigdf
