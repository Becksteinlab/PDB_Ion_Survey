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

def get_coord_num(ionname, atomname='O', bs=.1, mindistance=True, thres=.1):
    if not mindistance:
        df = pd.read_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins.csv')
    else:
        df = pd.read_csv(csvpath.abspath+'d-'+ionname+'-'+atomname+'-'+str(int(bs*100))+'pmbins-withmin.csv')
    locations = pd.read_csv(csvpath.abspath+'peaklocations-thres='+str(int(thres*10))+'.csv')
    heights = pd.read_csv(csvpath.abspath+'peakheights-thres='+str(int(thres*10))+'.csv')
    assert list(locations['ion/atom']).index("('{}', '{}')".format(ionname, atomname)) == list(heights['ion/atom']).index("('{}', '{}')".format(ionname, atomname))
    n = list(heights['ion/atom']).index("('{}', '{}')".format(ionname, atomname))
    m = max([heights['first max'][n], heights['second max'][n], heights['third max'][n]])
    if heights['first max'][n] == m and m > 0:
        rad = locations['first min'][n]
    elif heights['second max'][n] == m and m > 0:
        rad = locations['second min'][n]
    else:
        return None, None
    if not rad > 0:
        return None, None
    coordnum = 0
    for i in range(len(df)):
        if df['radius'][i] > rad:
            break
        else:
            mid = df['radius'][i]
            dm = df['radius'][1] - df['radius'][0]
            V = 4. / 3 * np.pi * ((mid+dm/2) ** 3 - (mid-dm/2) ** 3)
            d = df['density'][i]
            d *= V
            coordnum += d
    print(ionname, atomname, rad, coordnum)
    return rad, coordnum
