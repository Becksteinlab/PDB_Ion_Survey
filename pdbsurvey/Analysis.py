import mdsynthesis as mds
import MDAnalysis as mda
import coordination
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import peakutils

def make_sims(path):
    sims = mds.Tree('sims')
    pdbs = mds.Tree(path)
    pdbfiles = pdbs.glob('*.pdb')
    return pdbfiles

err_universe = []
def define_universe(pdbfile):
    s = mds.Sim(sims[os.path.splitext(pdbfile.name)[0] + '/'])
    with open(s[pdbfile.name].abspath, 'w') as f:
        f.write(pdbfile.read())
    try:
        s.universe = mda.Universe(s[pdbfile.name].abspath)
    except:
        with open('failures.out', 'a') as f:
            f.write(pdbfile.name + '\n')

def sim_labeling(sim):
    bundle.tags.add('pdb_ion_survey', 'pdbsurvey')

    f = open(sim.glob('*.pdb')[0].abspath, "r")
    searchlines = f.readlines()
    f.close()
    found = False
    for i, line in enumerate(searchlines):
        if ("REMARK" in line) and ("RESOLUTION." in line):
            try:
                things = line.split()
                resolutionindex = things.index('RESOLUTION.')
                value = things[resolutionindex + 1]
                value = float(value)
                found = True
                break
            except IndexError:
                continue
            except ValueError:
                continue
    if found:
        sim.categories['resolution'] = value    
    else:
        sim.categories['resolution'] = 'N/A'
    return value

    if sim.universe.select_atoms('resname HOH'):
        value = True
    else:
        value = False
    sim.categories['has_water'] = value

def closest_oxy_distance(bundle, ions, reses, cume=True):
    c = bundle
    for res in reses:
        if not cume:
            c = c[[r > (res - .5) for r in c.categories['resolution']]]
        c = c[[r <= res for r in c.categories['resolution']]]
        for ion in ions:
            z = c[c.tags[ion]]
            try:
                frames = []
                for d in z:
                    for iondata in d['coordination/' + ion.upper() + '/'].data:
                        frames.append(d['coordination/' + ion.upper() + '/'].data[iondata].sort('distance')[0:6].reset_index()['distance'])

                aoxy = []
                boxy = []
                coxy = []
                doxy = []
                eoxy = []
                foxy = []

                for i in range(len(frames)):
                    aoxy.append(frames[i][0])
                    boxy.append(frames[i][1])
                    coxy.append(frames[i][2])
                    doxy.append(frames[i][3])
                    eoxy.append(frames[i][4])
                    foxy.append(frames[i][5])

                bins = np.arange(0, 8, .2)

                ah, e = np.histogram(aoxy, bins=bins)
                am = .5 * (e[:-1] + e[1:])
                bh, e = np.histogram(boxy, bins=bins)
                bm = .5 * (e[:-1] + e[1:])
                ch, e = np.histogram(coxy, bins=bins)
                cm = .5 * (e[:-1] + e[1:])
                dh, e = np.histogram(doxy, bins=bins)
                dm = .5 * (e[:-1] + e[1:])
                eh, e = np.histogram(eoxy, bins=bins)
                em = .5 * (e[:-1] + e[1:])
                fh, e = np.histogram(foxy, bins=bins)
                fm = .5 * (e[:-1] + e[1:])

                fig = plt.figure(figsize = (4,3))
                ax = fig.add_subplot(1,1,1)
                ax.set_xlim(1, 6)

                ax.plot(am, ah, '#fe1e1e', label='1st oxy', lw=2)
                ax.plot(bm, bh, '#f7780f', label='2nd oxy', lw=2)
                ax.plot(cm, ch, '#39b927', label='3rd oxy', lw=2)
                ax.plot(dm, dh, '#1789fa', label='4th oxy', lw=2)
                ax.plot(em, eh, '#7f40ed', label='5th oxy', lw=2)
                ax.plot(fm, fh, '#ff2089', label='6th oxy', lw=2)

                ax.set_xlabel('Distance ($\AA$)')
                ax.set_ylabel('Frequency')
                ax.legend(fontsize='medium')

                if cume:
                    ax.figure.savefig(ion + "+ Distance of Oxygens Histogram (res from " + (res - .5) + " to " + res + ").pdf")
                else:
                    ax.figure.savefig(ion + "+ Distance of Oxygens Histogram (res to " + res + ").pdf")
            except:
                continue

def closest_oxy_distance(bundle, ions, reses, bins=200, cume=True):
    c = bundle
    for res in reses:
        if not cume:
            c = c[[r > (res - .5) for r in c.categories['resolution']]]
        c = c[[r <= res for r in c.categories['resolution']]]
        for ion in ions:
            z = c[c.tags[ion]]
            try:
                frames = []
                for d in z:
                    for iondata in d['coordination/' + ion.upper() + '/'].data:
                        frames.append(d['coordination/' + ion.upper() + '/'].data[iondata].sort('distance')[0:6].reset_index()['distance'])

                aoxy = []
                boxy = []
                coxy = []
                doxy = []
                eoxy = []
                foxy = []

                for i in range(len(frames)):
                    aoxy.append(frames[i][0])
                    boxy.append(frames[i][1])
                    coxy.append(frames[i][2])
                    doxy.append(frames[i][3])
                    eoxy.append(frames[i][4])
                    foxy.append(frames[i][5])

                ah, e = np.histogram(aoxy, bins=bins)
                am = .5 * (e[:-1] + e[1:])
                bh, e = np.histogram(boxy, bins=bins)
                bm = .5 * (e[:-1] + e[1:])
                ch, e = np.histogram(coxy, bins=bins)
                cm = .5 * (e[:-1] + e[1:])
                dh, e = np.histogram(doxy, bins=bins)
                dm = .5 * (e[:-1] + e[1:])
                eh, e = np.histogram(eoxy, bins=bins)
                em = .5 * (e[:-1] + e[1:])
                fh, e = np.histogram(foxy, bins=bins)
                fm = .5 * (e[:-1] + e[1:])

                fig = plt.figure(figsize = (4,3))
                ax = fig.add_subplot(1,1,1)
                ax.set_xlim(1, 6)

                ax.plot(am, ah, '#fe1e1e', label='1st oxy', lw=2)
                ax.plot(bm, bh, '#f7780f', label='2nd oxy', lw=2)
                ax.plot(cm, ch, '#39b927', label='3rd oxy', lw=2)
                ax.plot(dm, dh, '#1789fa', label='4th oxy', lw=2)
                ax.plot(em, eh, '#7f40ed', label='5th oxy', lw=2)
                ax.plot(fm, fh, '#ff2089', label='6th oxy', lw=2)

                ax.set_xlabel('Distance ($\AA$)')
                ax.set_ylabel('Frequency')
                ax.legend(fontsize='medium')

                if cume:
                    ax.figure.savefig(ion + "+ Distance of Oxygens Histogram (res from " + (res - .5) + " to " + res + ").pdf")
                else:
                    ax.figure.savefig(ion + "+ Distance of Oxygens Histogram (res to " + res + ").pdf")
            except:
                continue
