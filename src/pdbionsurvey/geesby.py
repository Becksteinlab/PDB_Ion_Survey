from __future__ import division
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
# import pdbionsurvey.analysis
from os import path as pth
import os
import shutil
import peakutils
from glob import glob
# %matplotlib inline

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                       AutoMinorLocator)
import seaborn as sns; sns.set()
sns.set_palette(sns.color_palette('husl'))
import mmtf
import warnings

plt.style.use('ggplot')
sns.set_style('ticks')

def make_gees_bytech_df(ionname, atomname='O', bs=.1, methods=None, label='total', bundle=b, path=csvpath):
    """
    :Arguments:
    """
    gdfs = {}
    if not methods:
        methods = set(bundle.categories['method'])
    bundle = bundle[bundle.tags[ionname]]
    raddone = False
    for tech in methods:
        try:
            tempbundle = bundle[bundle.tags[tech]]
            gdf = pdbionsurvey.coordination.gee(tempbundle, ionname, atomname=atomname, binsize=bs)
            gdf['density'] = gdf['density']/bulkdensity[atomname]
            gdfs[tech] = gdf['density']
            if not raddone:
                gdfs['radius'] = gdf['radius']
                raddone = True
        except ValueError:
            pass
        print(tech, 'done')
    df = pd.DataFrame(gdfs)
    df.to_csv(path.abspath+'g-bytech-'+ionname+'-'+atomname+'('+label+').csv')
    return df

def make_gees_bytech(ionname, atomname='O', bs=.1, maxdistance=15, mindistance=.5, methods=None, label='total', bundle=b, readpath=csvpath, path=impath):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    if not methods:
        methods = set(bundle.categories['method'])
    gdf = pd.read_csv(path.abspath+'g-bytech-'+ionname+'-'+atomname+'('+label+').csv')
    gdf = gdf[gdf['radius'] < maxdistance]
    ys = []

    for tech in methods:
        try:
            gdf[tech][gdf['radius'] < mindistance] = 0
            ys.append(max(gdf[tech]))
            if atomname == 'C':
                ax.plot(gdf['radius'], gdf[tech], label=tech, alpha=.5, linewidth=2)
            else:
                ax.plot(gdf['radius'], gdf[tech], label=tech, linewidth=2)
        except KeyError:
            continue
    ts = getbins(ax.get_xlim()[1])
    ax.xaxis.set_major_locator(MultipleLocator(5*ts))
    ax.xaxis.set_minor_locator(MultipleLocator(ts))
    sns.despine(offset=10, ax=ax)
    ax.set_xlabel(r'distance ($\mathrm{\AA}$)')
    ax.set_ylabel(r'$g(r)$')
    yts = getbins(max(ys))
    ax.yaxis.set_major_locator(MultipleLocator(5*yts))
    ax.yaxis.set_minor_locator(MultipleLocator(yts))
    ax.set_xlim(0, maxdistance)
    ax.set_ylim(0, max(ys))
    ax.plot(ax.get_xlim(), [1, 1], color=(0,0,0), ls='dotted', alpha=.5)
    ax.legend()
    ax.figure.savefig(path.abspath+'g-bytech-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-('+label+').png')
    ax.figure.savefig(path.abspath+'g-bytech-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-('+label+').pdf')
    return ax

x = [['NA', 'K', 'LI'], ['CL', 'F', 'BR'], ['ZN', 'MG', 'CA'], ['MN', 'CU', 'BA'], ['CD', 'HG', 'TL']]

def make_gees_bylig(ionnames, maxdistance=5, bs=.1, ts=.1):
    for atomname in ['O', 'N', 'S', 'C']:
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        fig.set_tight_layout(True)
        fig1 = plt.figure(figsize=(4,3))
        ax1 = fig1.add_subplot(111)
        fig1.set_tight_layout(True)
        ys = []
        for ionname in ionnames:
            print('started g '+ionname+' with '+atomname)
            gdf = pd.read_csv('d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-withmin.csv')
            gdf = gdf[gdf['radius'] < maxdistance]
            gdf['density'] = [gdf['density'][i] if gdf['radius'][i]<.5 else 0 for i in range(len(gdf['density']))]
            ax.plot(gdf['radius'], gdf['density'], label=propernames[ionname], linewidth=2)
            y = gdf['density']/bulkdensity[atomname]
            ax1.plot(gdf['radius'], y, label=propernames[ionname], linewidth=2)
            ys.append(max(y))
        ax.set_xlabel(r'distance ($\mathrm{\AA}$)')
        ax.set_ylabel(r'density ($\mathrm{\AA}^{-3}$)')
        ax.xaxis.set_major_locator(MultipleLocator(5*ts))
        ax.xaxis.set_minor_locator(MultipleLocator(ts))
        sns.despine(offset=10, ax=ax)
        ax.legend()
        charge = 'random'
        if all([pdbionsurvey.coordination.get_charge(ionname) == pdbionsurvey.coordination.get_charge(ionnames[0]) for ionname in ionnames]):
            c = pdbionsurvey.coordination.get_charge(ionnames[0])
            if c < 0:
                charge = 'minus'+str(c)
            else:
                charge = 'plus'+str(c)
        ax.figure.savefig(impath.abspath+'absolutedensity-'+charge+'_'.join(ionnames)+'-'+atomname+'-'+str(int(bs*100))+'pmbins'+str(maxdistance)+'.png')
        ax.figure.savefig(impath.abspath+'absolutedensity-'+charge+'_'.join(ionnames)+'-'+atomname+'-'+str(int(bs*100))+'pmbins'+str(maxdistance)+'.pdf')
        ax1.set_xlabel(r'distance ($\mathrm{\AA}$)')
        ax1.set_ylabel(r'$g(r)$')
        ax1.xaxis.set_major_locator(MultipleLocator(1))
        ax1.xaxis.set_minor_locator(MultipleLocator(ts))
        yts = getbins(max(ys))
        ax1.yaxis.set_major_locator(MultipleLocator(yts))
        ax1.yaxis.set_minor_locator(MultipleLocator(5*yts))
        ax1.set_xlim(0, maxdistance)
        ax1.set_ylim(0, max(ys))
        ax1.plot(gdf['radius'], np.array([1 for i in range(len(gdf['radius']))]), color=(0,0,0), ls='dotted', alpha=.5)
        sns.despine(offset=10, ax=ax1)
        ax1.legend()
        ax1.figure.savefig(impath.abspath+'g-'+charge+'-'+'_'.join(ionnames)+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.png')
        ax1.figure.savefig(impath.abspath+'g-'+charge+'-'+'_'.join(ionnames)+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.pdf')

def make_gees_byion(ionname, atomnames=['O', 'N', 'S', 'C'], maxdistance=5, bs=.1, ts=.1):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    ys = []
    for atomname in atomnames:
        print('started g '+ionname+' with '+atomname)
        gdf = pd.read_csv(path.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-withmin.csv')
        gdf = gdf[gdf['radius'] < maxdistance]
        gdf['density'] = [gdf['density'][i] if gdf['radius'][i]<.5 else 0 for i in range(len(gdf['density']))]
        y = gdf['density']/bulkdensity[atomname]
        if atomname == 'C':
            ax.plot(gdf['radius'], gdf['density'], label=atomname, alpha=.5, linewidth=2)
            ax1.plot(gdf['radius'], y, label=atomname, alpha=.5, linewidth=2)
        else:
            ax.plot(gdf['radius'], gdf['density'], label=atomname, linewidth=2)
            ax1.plot(gdf['radius'], y, label=atomname, linewidth=2)
        ys.append(max(y))
    ax.set_xlabel(r'distance ($\mathrm{\AA}$)')
    ax.set_ylabel(r'density ($\mathrm{\AA}^{-3}$)')
    ts = getbins(max(gdf['radius']))
    ax.xaxis.set_major_locator(MultipleLocator(5*ts))
    ax.xaxis.set_minor_locator(MultipleLocator(ts))
    sns.despine(offset=10, ax=ax)
    ax.legend()
    ax.figure.savefig(impath.abspath+'absolutedensity-'+ionname+'-'+'_'.join(atomnames)+'-'+str(int(bs*100))+'pmbins'+str(maxdistance)+'.png')
    ax.figure.savefig(impath.abspath+'absolutedensity-'+ionname+'-'+'_'.join(atomnames)+'-'+str(int(bs*100))+'pmbins'+str(maxdistance)+'.pdf')
    ax1.set_xlabel(r'distance ($\mathrm{\AA}$)')
    ax1.set_ylabel(r'$g(r)$')
    ts = getbins(max(gdf['radius']))
    ax1.xaxis.set_major_locator(MultipleLocator(5*ts))
    ax1.xaxis.set_minor_locator(MultipleLocator(ts))
    yts = getbins(max(ys))
    ax1.yaxis.set_major_locator(MultipleLocator(5*yts))
    ax1.yaxis.set_minor_locator(MultipleLocator(yts))
    ax1.set_xlim(0, maxdistance)
    ax1.set_ylim(0, max(ys))
    ax1.plot(gdf['radius'], np.array([1 for i in range(len(gdf['radius']))]), color=(0,0,0), ls='dotted', alpha=.5)
    sns.despine(offset=10, ax=ax1)
    ax1.legend()
    ax1.figure.savefig(impath.abspath+'g-'+ionname+'-'+'_'.join(atomnames)+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.png')
    ax1.figure.savefig(impath.abspath+'g-'+ionname+'-'+'_'.join(atomnames)+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.pdf')
    return ax, ax1

def make_gees_byiongroups(ionname, binsizes=[.1], maxes=[5, 10, 15]):
    for md in maxes:
        for bs in binsizes:
            plot_together_byion(ionname, atomnames=['O', 'N'], maxdistance=md, bs=bs, ts=md/25)
            plot_together_byion(ionname, atomnames=['O', 'N', 'S'], maxdistance=md, bs=bs, ts=md/25)
            plot_together_byion(ionname, atomnames=['O', 'N', 'S', 'C'], maxdistance=md, bs=bs, ts=md/25)

def make_gees_byres_df(ionname, atomname='O', reses=RESES, bundle=b, path=csvpath):
    bs=.1
    naexists = True
    bundle = bundle[bundle.tags[ionname]]
    resna = bundle[[r=='N/A' for r in bundle.categories['resolution']]]
    bund = bundle - resna
    reses = list(np.sort(reses)).reverse()
    gdfs = {}
    try:
        gdfany = pdbionsurvey.coordination.gee(bundle, ionname, atomname=atomname, binsize=bs)
        gdfany['density'] = gdfany['density']/bulkdensity[atomname]
        gdfs['any'] = gdfany['density']
    except ValueError:
        pass
    print('any done')
    try:
        gdfna = pdbionsurvey.coordination.gee(resna, ionname, atomname=atomname, binsize=bs)
        gdfna['density'] = gdfna['density']/bulkdensity[atomname]
        gdfs['unknown'] = gdfna['density']
    except ValueError:
        naexists = False
        pass
    print('na done')
    try:
        gdfnotna = pdbionsurvey.coordination.gee(bund, ionname, atomname=atomname, binsize=bs)
        gdfnotna['density'] = gdfnotna['density']/bulkdensity[atomname]
        gdfs['known'] = gdfnotna['density']
    except ValueError:
        pass
    print('not na done')
    for res in reses:
        try:
            gdf = pdbionsurvey.coordination.gee(bundle, ionname, atomname=atomname, binsize=bs)
            gdf['density'] = gdf['density']/bulkdensity[atomname]
            gdfs['res<='+str(res)] = gdf['density']
        except ValueError:
            pass
    print('reses done')
    gdfs['radius'] = gdfany['radius']
    df = pd.DataFrame(gdfs)
    df.to_csv(path.abspath+'g-byres-'+ionname+'-'+atomname+'-(any_known_unknown'+'_'.join(str(reses))+').csv')
    return df, naexists

def make_gees_byres(ionname, atomname='O', reses = RESES, maxdistance=5, mindistance=.5, bs=.1, bundle=b, path=impath):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    gdf, naexists = make_gees_byres_df(ionname, atomname=atomname, reses=reses, bundle=bundle)
    gdf = gdf[gdf['radius'] < maxdistance]
    ys = []
    if naexists:
        reslabels = ['any', 'unknown', 'known'] + ['res<='+str(res) for res in reses]
    else:
        reslabels = ['any'] + ['res<='+str(res) for res in reses]
    for reslabel in reslabels:
        try:
            gdf[reslabel][gdf['radius'] < mindistance] = 0
            ys.append(max(gdf[reslabel]))
            if atomname == 'C':
                ax.plot(gdf['radius'], gdf[reslabel], label=reslabel, alpha=.5, linewidth=2)
            else:
                ax.plot(gdf['radius'], gdf[reslabel], label=reslabel, linewidth=2)
        except KeyError:
            continue
    ts = getbins(ax.get_xlim()[1])
    ax.xaxis.set_major_locator(MultipleLocator(5*ts))
    ax.xaxis.set_minor_locator(MultipleLocator(ts))
    sns.despine(offset=10, ax=ax)
    ax.set_xlabel(r'distance ($\mathrm{\AA}$)')
    ax.set_ylabel(r'$g(r)$')
    yts = getbins(max(ys))
    ax.yaxis.set_major_locator(MultipleLocator(5*yts))
    ax.yaxis.set_minor_locator(MultipleLocator(yts))
    ax.set_xlim(0, maxdistance)
    ax.set_ylim(0, max(ys))
    ax.plot(ax.get_xlim(), [1, 1], color=(0,0,0), ls='dotted', alpha=.5)
    sns.despine(offset=10, ax=ax1)
    ax.legend()
    ax.figure.savefig(impath.abspath+'g-byres-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-(any_known_unknown'+'_'.join(str(reses))+').png')
    ax.figure.savefig(impath.abspath+'g-byres-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-(any_known_unknown'+'_'.join(str(reses))+').pdf')
