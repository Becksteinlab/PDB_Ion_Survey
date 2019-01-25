# PDB Ion Survey 
# Copyright (c) 2016 Kacey Clark
# Published under the GPL v3
# https://github.com/Becksteinlab/PDB_Ion_Survey/

'''
Functions for analyzing coordination functions
'''

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

def kolmogorov_smirnov(n1, n2, ylabel='count', cumulative=False):
    '''Gives the kolmogorov-smirnov statistic between two cumulative normalized discrete functions.
    :Arguments:
        *n1*
            pandas.DataFrame first n(r), order does not matter
        *n2*
            pandas.DataFrame second n(r), order does not matter
        *ylabel*
            String label of axis of interest in n1 and n2; default='count'
        *cumulative*
            Boolean true if n1 and n2 are already cumulative; default=False
    :Returns:
        *maxdiff*
            Float value of kolmogorov-smirnov statistic
    '''
    assert(np.equal(np.array(g1['radius']), np.array(g2['radius'])).all())
    if not cumulative:
        n1d = np.cumsum(np.array(n1[ylabel]))
        n2d = np.cumsum(np.array(n2[ylabel]))
    else:
        n1d = np.array(n1[ylabel])
        n2d = np.array(n2[ylabel])
    n1d = n1['count']/max(n1['count'])
    n2d = n2['count']/max(n2['count'])
    diff = abs(n2d - n1d)
    maxdiff = max(diff)
    return maxdiff

def cramer_vonmises(n1, n2, ylabel='count', cumulative=False):
    '''Gives the smirnov-cramer-vonmises statistic between two cumulative normalized discrete functions.
    :Arguments:
        *n1*
            pandas.DataFrame first n(r), order does not matter
        *n2*
            pandas.DataFrame second n(r), order does not matter
        *ylabel*
            String label of axis of interest in n1 and n2; default='count'
        *cumulative*
            Boolean true if n1 and n2 are already cumulative; default=False
    :Returns:
        *wsquared*
            Float value of cramer-von_mises statistic
    '''
    assert(np.equal(np.array(n1['radius']), np.array(n2['radius'])).all())
    if not cumulative:
        n1d = np.cumsum(np.array(n1[ylabel]))
        n2d = np.cumsum(np.array(n2[ylabel]))
    else:
        n1d = np.array(n1[ylabel])
        n2d = np.array(n2[ylabel])
    diff = abs(n2d - n1d)
    squarediff = np.square(diff)
    parts = squarediff*(g1d-np.insert(g1d[:-1], 0, [0]))
    wsquared = np.sum(parts)
    return wsquared

def c(alpha):
    '''Gives the kolmogorov-smirnov statistic between two cumulative normalized discrete functions.
    :Arguments:
        *alpha*
            Float confidence
    :Returns:
        *c*
            Float c
    '''
    c = np.sqrt(-1/2*np.log(alpha))
    return c

def kolmogorov_smirnov_test(g1, g2, ylabel='count', cumulative=False, alpha=.001):
    '''Gives the result of the kolmogorov-smirnov test of two cumulative normalized discrete functions.
    :Arguments:
        *g1*
            pandas.DataFrame first n(r), order does not matter
        *g2*
            pandas.DataFrame second n(r), order does not matter
        *ylabel*
            String label of axis of interest in n1 and n2; default='count'
        *cumulative*
            Boolean true if n1 and n2 are already cumulative; default=False
        *alpha*
            Float confidence; default=.1
    :Returns:
        *res*
            Boolean true if g1, g2 distinct with alpha confidence
    '''
    ks = kolmogorov_smirnov(g1, g2, ylabel='count', cumulative=False)
    n1 = len(g1)
    n2 = len(g2)
    calpha = c(alpha)
    std = calpha * np.sqrt((n1+n2)/(n1*n2))
    if ks > std:
        print('Distinct with ', alpha, 'confidence.')
        res= True
    else:
        print('Not distinct with', alpha, 'confidence.')
        res= False
    return res

def ks_test_comparison(ionnames=IONNAMES, alpha=1):
    '''Gives the result of the kolmogorov-smirnov test of g(r)s of O, N, S around many ions.
    :Arguments:
        *ionnames*
            List ionnames to compare; default=IONNAMES
        *alpha*
            Float confidence; default=.1
    :Returns:
        *bigdf*
            pd.DataFrame dataframe of kolmogorov-smirnov test results for different ion/atom combinations
    '''
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
