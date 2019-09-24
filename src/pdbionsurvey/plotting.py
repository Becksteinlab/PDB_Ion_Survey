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

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                       AutoMinorLocator)
import seaborn as sns; sns.set()
sns.set_palette(sns.color_palette('husl'))
import mmtf

plt.style.use('ggplot')
sns.set_style('ticks')

#table 1 in Peters paper
newcations = ['NA', 'MG', 'K', 'CA', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'PD', 'AG', 'CD', 'IR', 'PT', 'AU', 'HG']

othercations = ['LA', 'PB', 'TL', 'LI', 'BA', 'RB', 'CS', 'SR']

anions = ['CL', 'F', 'BR']

omit = ['V', 'IOD']

IONNAMES = list(set(newcations).union(set(othercations), set(anions)))

propernames = {'AG':r'Ag$^{+}$', 'AU': r'Au$^{+}$', 'BA': r'Ba$^{2+}$',
               'BR': r'Br$^{-}$', 'CA': r'Ca$^{2+}$', 'CD': r'Cd$^{2+}$',
               'CL': r'Cl$^{-}$', 'CO': r'Co$^{2+}$', 'CR': r'Cr$^{3+}$',
               'CS': r'Cs$^{+}$', 'CU': r'Cu$^{2+}$', 'FE':'Fe$^{3+}$',
               'F': r'F$^{-}$', 'HG': r'Hg$^{2+}$', 'IR': r'Ir$^{4+}$',
               'IOD': r'I$^{-}$', 'K': r'K$^{+}$', 'LA': 'La$^{3+}$',
               'LI': r'Li$^{+}$', 'MG': r'Mg$^{2+}$', 'MN': r'Mn$^{2+}$',
               'NA':r'Na$^{+}$', 'NI': r'Ni$^{2+}$', 'PB': r'Pb$^{2+}$',
               'PD': r'Pd$^{2+}$', 'PT': r'Pt$^{2+}$', 'RB': r'Rb$^{+}$',
               'SR': r'Sr$^{2+}$', 'TL': r'Tl$^{+}$', 'V': r'V$^{3+}$',
               'ZN':r'Zn$^{2+}$'}

ATOMNAMES = ['O', 'N', 'S', 'C']

sel = {'O': 'name O* and not name OS', 'N': 'name N* and not name NA and not name NE', 'S': 'name S* and not name SB', 'C': 'name C* and not name CS and not resname CA'}

obulk = 0.00947
nbulk = 0.00609
sbulk = 0.000196
cbulk = 0.0234

bulkdensity = {'O': obulk, 'N': nbulk, 'S': sbulk, 'C': cbulk}

# directory for the pdbs
PATH = dtr.Tree('/nfs/homes3/kreidy/Projects/pdbsurvey/')
csvpath = PATH['csvs/']
impath = PATH['images/']

# directory for the sims
# it is important that this directory is empty,
# otherwise, mislabeling may happen later
protdir = PATH['proteins/']

# directory for DataFrames
# datapath = path['dataframes']

RESES = [1, 1.4, 2]

b = dtr.discover(protdir)

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

def make_dees(ionname, atomnames=ATOMNAMES, bs=.1, maxdistance=15):
    for atomname in atomnames:
        print('started g '+ionname+' with '+atomname)
        gdf = pdbionsurvey.coordination.gee(b, ionname, atomname=atomname, binsize=bs)
        gdf = gdf[gdf['radius'] < maxdistance]
        print('made g '+ionname+' with '+atomname)
        gdf.to_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
        print('saved g '+ionname+' with '+atomname)
#         fig = plt.figure(figsize=(4,3))
#         ax = fig.add_subplot(111)
#         fig.set_tight_layout(True)
#         ax.plot(gdf['radius'], gdf['density'], linewidth=2)
#         print('plotted g '+ionname+' with '+atomname)
#         ax.set_xlabel(r'distance ($\mathrm{\AA}$)')
#         ax.set_ylabel(r'density ($\mathrm{\AA} ^{-3}$)')
#         ax.set_xticks(np.arange(0, maxdistance+ts, ts))
#     #         ax.set_title('Average density RDF of '+atomname+' around '+ionname)
#         sns.despine(offset=10, ax=ax)
#         ax.figure.savefig(path.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.png')
#         ax.figure.savefig(path.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'.pdf')
#         print('finished g '+ionname+' with '+atomname)

def make_gees(ionname, atomname='O', maxdistance=5, bs=.1, bundle=b, path=csvpath):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    print('started g '+ionname+' with '+atomname)
    mindistance = .5
    gdf = pd.read_csv(path.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
    gdf['density'][gdf['radius'] < mindistance] = 0
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
    print('finished g '+ionname+' with '+atomname)

methods = set(b.categories['method'])

METHODS = list(methods)

def make_gees_bytech_df(ionname, atomname='O', bs=.1, methods=METHODS, label='total', bundle=b, path=csvpath):
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
    gdf = pd.read_csv(readpath.abspath+'g-bytech-'+ionname+'-'+atomname+'('+label+').csv')
    gdf = gdf[gdf['radius'] < maxdistance]
    ys = []

    techcolors = {'FIBER DIFFRACTION': '#93F388', 'POWDER DIFFRACTION': '#F576B9', 'X-RAY DIFFRACTION': '#0208B2',
                   'SOLID-STATE NMR': '#DE28E7', 'ELECTRON CRYSTALLOGRAPHY': '#16AAB9', 'SOLUTION NMR': '#F12A0F',
                   'NEUTRON DIFFRACTION': '#E3F10F', 'multiple': '#6F16B9', 'ELECTRON MICROSCOPY': '#0BA032'}

    for tech in methods:
        try:
            gdf[tech][gdf['radius'] < mindistance] = 0
            ys.append(max(gdf[tech]))
            ax.plot(gdf['radius'], gdf[tech], label=tech, linewidth=2, color=techcolors[tech])
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
    ax.figure.savefig(path.abspath+'g-bytech-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-('+label+').png')
    ax.figure.savefig(path.abspath+'g-bytech-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-('+label+').pdf')
    return ax

def make_gees_bytech_df_smallbunds(ionname, atomname='O', bs=.1, methods=None, label='total', bundle=b, path=csvpath):
    if not methods:
        methods = set(bundle.categories['method'])
    bund = bundle
    for i in range(0, len(bund), 1000):
        gdfs = {}
        bundle = bund[i, i+1000]
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

        df.to_csv(path.abspath+'g(small)-bytech-'+ionname+'-'+atomname+'('+str(i)+').csv')
    return df

def make_gees_bytech_smallbunds(ionname, atomname='O', bs=.1, maxdistance=15, mindistance=.5, methods=None, label='total', bundle=b, readpath=csvpath, path=impath):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    df = {}
    if not methods:
        methods = set(bundle.categories['method'])
    for i in range(0, len(bundle), 1000):
        bund = bundle[i, i+1000]
        gdf = pd.read_csv(readpath.abspath+'g(small)-bytech-'+ionname+'-'+atomname+'('+str(i)+').csv')
        if i == 0:
            df['radius'] = gdf['radius']
            for tech in methods:
                df[tech] = len(bund)/len(bundle)*gdf[tech]
        else:
            for tech in methods:
                df[tech] += len(bund)/len(bundle)*gdf[tech]
    gdf = df[df['radius'] < maxdistance]
    ys = []

    techcolors = {'FIBER DIFFRACTION': '#93F388', 'POWDER DIFFRACTION': '#F576B9', 'X-RAY DIFFRACTION': '#0208B2',
                   'SOLID-STATE NMR': '#DE28E7', 'ELECTRON CRYSTALLOGRAPHY': '#16AAB9', 'SOLUTION NMR': '#F12A0F',
                   'NEUTRON DIFFRACTION': '#E3F10F', 'multiple': '#6F16B9', 'ELECTRON MICROSCOPY': '#0BA032'}

    for tech in methods:
        try:
            gdf[tech][gdf['radius'] < mindistance] = 0
            ys.append(max(gdf[tech]))
            ax.plot(gdf['radius'], gdf[tech], label=tech, linewidth=2, color=techcolors[tech])
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
    ax.figure.savefig(path.abspath+'g-bytech-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-('+label+').png')
    ax.figure.savefig(path.abspath+'g-bytech-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-('+label+').pdf')
    return ax

x = [['NA', 'K', 'LI'], ['CL', 'F', 'BR'], ['ZN', 'MG', 'CA'], ['MN', 'CU', 'BA'], ['CD', 'HG', 'TL']]

def plot_together_bylig(ionnames, maxdistance=5, bs=.1):
    for atomname in ['O', 'N', 'S', 'C']:
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        fig.set_tight_layout(True)
        fig1 = plt.figure(figsize=(4,3))
        ax1 = fig1.add_subplot(111)
        fig1.set_tight_layout(True)
        ys = []
        mindistance=.5
        for ionname in ionnames:
            print('started g '+ionname+' with '+atomname)
            gdf = pd.read_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
            gdf = gdf[gdf['radius'] < maxdistance]
            gdf['density'][gdf['radius'] < mindistance] = 0
            ax.plot(gdf['radius'], gdf['density'], label=propernames[ionname], linewidth=2)
            y = gdf['density']/bulkdensity[atomname]
            ax1.plot(gdf['radius'], y, label=propernames[ionname], linewidth=2)
            ys.append(max(y))
        ax.set_xlabel(r'distance ($\mathrm{\AA}$)')
        ax.set_ylabel(r'density ($\mathrm{\AA}^{-3}$)')
        ts = getbins(ax.get_xlim()[1])
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

def plot_together_byion(ionname, atomnames=['O', 'N', 'S', 'C'], maxdistance=5, bs=.1, path=csvpath):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    ys = []
    mindistance = .5
    for atomname in atomnames:
        print('started g '+ionname+' with '+atomname)
        gdf = pd.read_csv(path.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
        gdf = gdf[gdf['radius'] < maxdistance]
        gdf['density'][gdf['radius'] < mindistance] = 0
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

def plot_together_byion_onaxes(ionname, atomnames=['O', 'N', 'S', 'C'], maxdistance=5, bs=.1, path=csvpath):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    fig.set_tight_layout(True)
    fig1 = plt.figure(figsize=(4,3))
    ax1 = fig1.add_subplot(111)
    fig1.set_tight_layout(True)
    ys = []
    mindistance = .5
    for atomname in atomnames:
        print('started g '+ionname+' with '+atomname)
        gdf = pd.read_csv(path.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
        gdf = gdf[gdf['radius'] < maxdistance]
        gdf['density'][gdf['radius'] < mindistance] = 0
        y = gdf['density']/bulkdensity[atomname]
        if atomname == 'C':
            ax.plot(gdf['radius'], gdf['density'], label=atomname, alpha=.5, linewidth=2)
            ax1.plot(gdf['radius'], y, label=atomname, alpha=.5, linewidth=2)
        else:
            ax.plot(gdf['radius'], gdf['density'], label=atomname, linewidth=2)
            ax1.plot(gdf['radius'], y, label=atomname, linewidth=2)
        if atomname == 'O' or atomname == 'N':
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

    if reses is not None:
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
        if reses is not None:
            reslabels = ['any', 'unknown', 'known'] + ['res<='+str(res) for res in reses]
        else:
            reslabels = ['any', 'unknown', 'known']
    else:
        if reses is not None:
            reslabels = ['any'] + ['res<='+str(res) for res in reses]
        else:
            reslabels = ['any']
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
    ax.figure.savefig(path.abspath+'g-byres-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-(any_known_unknown'+'_'.join(str(reses))+').png')
    ax.figure.savefig(impath.abspath+'g-byres-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-'+str(maxdistance)+'-(any_known_unknown'+'_'.join(str(reses))+').pdf')
