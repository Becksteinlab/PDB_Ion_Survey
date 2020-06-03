
m __future__ import division
import MDAnalysis as mda
import datreant as dtr
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pdbionsurvey.coordination
import json
import seaborn as sb
import scipy
import mmtf
import pdbionsurvey.collection
from matplotlib.ticker import MaxNLocator
# import pdbionsurvey.analysis
from os import path as pth
import os
import shutil
from glob import glob

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                       AutoMinorLocator)
import mmtf

def individualdee(prot, ionname='NA', atomname='O', bs=.1, mindistance=True, maxdistance=15, ts=1):
    u = mda.Universe(prot[prot.name+'.pdb'].abspath)
    ions = u.select_atoms('name '+ionname+' and resname '+ionname)
    for ion in ions:
        frames = []
        for csv in prot.glob('coordination/'+ionname.upper()+'/'+atomname+'/*.csv'):
            df = pd.read_csv(csv.abspath)
            frames.append(df)
        dataframe = pd.concat(frames)

        h, e = np.histogram(dataframe['distance'], bins=np.arange(0, max(dataframe['distance']), bs))
        m = .5 * (e[:-1] + e[1:])
        V = 4. / 3 * np.pi * (e[1:] ** 3 - e[:-1] ** 3)
        density = h / V

        gdf = pd.DataFrame({'radius': m, 'density': density}, columns=['radius', 'density'])

        gdf = gdf[gdf['radius'] < maxdistance]
        if not mindistance:
            gdf.to_csv(csvpath.abspath+'individual/d-'+prot.name+'-'+ionname.upper()+'-'+str(ion.resnum)+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
        else:
            mindistance = .5
            gdf['density'] = [gdf['density'][i] if gdf['radius'][i]>.5 else 0 for i in range(len(gdf['density']))]
            gdf.to_csv(csvpath.abspath+'individual/d-'+prot.name+'-'+ionname.upper()+'-'+str(ion.resnum)+'-'+atomname+'-'+str(int(bs*100))+'pmbins-withmin.csv')

def make_dees(ionname, atomnames=ATOMNAMES, bs=.1, mindistance=True, maxdistance=15, ts=1):
    for atomname in atomnames:
        print('started d '+ionname+' with '+atomname)
        gdf = pdbionsurvey.coordination.gee(b, ionname, atomname=atomname, binsize=bs)
        gdf = gdf[gdf['radius'] < maxdistance]
        print('made d '+ionname+' with '+atomname)
        if not mindistance:
            gdf.to_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
        else:
            mindistance = .5
            gdf['density'] = [gdf['density'][i] if gdf['radius'][i]>mindistance else 0 for i in range(len(gdf['density']))]
            gdf.to_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-withmin.csv')
        print('saved d '+ionname+' with '+atomname)

def getbins(num):
    ts = int(num/10)
    if ts == 0:
        ts = num/10
    if ts < 1 and ts >= .5:
        ts = .5
    elif ts < .5 and ts >= .2:
        ts = .2
    elif ts < .2:
        ts = .1
    return ts

def make_gees(ionname, atomname='O', maxdistance=15, bs=.1, bundle=b, path=csvpath):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    print('started g '+ionname+' with '+atomname)
    gdf = pd.read_csv(path.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-withmin.csv')
    gdf['density'] = [gdf['density'][i] if gdf['radius'][i]<.5 else 0 for i in range(len(gdf['density']))]
    gdf = gdf[gdf['radius'] < maxdistance]
    ax.plot(gdf['radius'], gdf['density'], label=propernames[ionname], linewidth=2)

    yts = getbins(max(gdf['density']))
    ts = getbins(max(gdf['radius']))

    ax.set_xlabel(r'distance ($\mathrm{\AA}$)')
    ax.set_ylabel(r'density ($\mathrm{\AA}^{-3}$)')
    ax.xaxis.set_major_locator(MultipleLocator(5*ts))
    ax.xaxis.set_minor_locator(MultipleLocator(ts))
    ax1.yaxis.set_major_locator(MultipleLocator(.005))
    ax1.yaxis.set_minor_locator(MultipleLocator(yts*.005))
    sns.despine(offset=10, ax=ax)
    ax.legend()

    ax.figure.savefig(impath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.png')
    ax.figure.savefig(impath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.pdf')

    y = gdf['density']/bulkdensity[atomname]
    ax1.plot(gdf['radius'], y, label=propernames[ionname], linewidth=2)

    yts = getbins(max(y))
    ts = getbins(max(gdf['radius']))

    ax1.set_xlabel(r'distance ($\mathrm{\AA}$)')
    ax1.set_ylabel(r'$g(r)$')
    ax1.xaxis.set_major_locator(MultipleLocator(5*ts))
    ax1.xaxis.set_minor_locator(MultipleLocator(ts))
    ax1.yaxis.set_major_locator(MultipleLocator(5*yts))
    ax1.yaxis.set_minor_locator(MultipleLocator(yts))
    ax1.set_xlim(0, maxdistance)
    ax1.set_ylim(0, max(y))
    sns.despine(offset=10, ax=ax1)
    ax1.legend()

    ax1.plot(gdf['radius'], np.array([1 for i in range(len(gdf['radius']))]), color=(0,0,0), ls='dotted', alpha=.5)

    ax1.figure.savefig(impath.abspath+'g-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.png')
    ax1.figure.savefig(impath.abspath+'g-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.pdf')
    df = pd.DataFrame({'radius': gdf['radius'], 'density': y}, columns=['radius', 'density'])
    df.to_csv(csvpath.abspath+'g-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.csv')
    print('finished g '+ionname+' with '+atomname)

IONNAME = 'CL'
ATOMNAME = 'O'
def individualen(prot, ionname=IONNAME, atomname=ATOMNAME, bs=.1, mindistance=True, ts=1):
    shellsize = shells['first min'][(IONNAME, ATOMNAME)]
    maxdistance = shellsize
    u = mda.Universe(prot[prot.name+'.pdb'].abspath)
    coordnums = []
    for csv in prot.glob('coordination/'+ionname.upper()+'/'+atomname+'/*.csv'):
        df = pd.read_csv(csv.abspath)
        try:
            gdf = df[df['distance'] < maxdistance]
            if mindistance:
                gdf = gdf[gdf['distance'] > .5]
            coordnum = len(gdf['distance'])
        except TypeError:
            coordnum = 'N/A'
            continue
        coordnums.append(coordnum)
        atomnum = csv.name[:-4]
        prot.categories['coordnum_'+atomnum] = coordnum
    prot.categories[IONNAME+'_'+ATOMNAME+'_coordnums'] = str(coordnums)
    return coordnums
