## service module to explore data and plot figures.
# pylint: disable=no-self-argument, no-member

from covid19 import get_curves_from, get_subjects_from, get_subject_number_from, get_vars_from, get_stats_from, show_, AUC_of
from core import value as dataset
from core import struct
from os import listdir

import numpy as np
from scipy import stats

from matplotlib import pyplot as lab
from matplotlib import colors

key = {
    'goHome' : 'go',
    'MVbio'  : 'mv',
    'ICUbio' : 'ic',
    'dcr'    : 'd',
    'lab'    : 'l',
    'CoV2'   : 'c',
}

def count_subjects_from(data):
    base = get_subjects_from(data)
    if type(base) is tuple: return list(set(base[0]+base[1]))
    return list(set(base))

#use this from an ipython console (easier way) after navigating to results folder
#the returned structure is explorable using structure methods and properties, see README.txt in results
def load(key, ftype='.res'):
    ''' load all results and return them as a structure: e.g., R = load(key) '''
    R = struct()
    for file in listdir():
        if file.endswith(ftype) and not 'comp' in file:
            name = file[:-len(ftype)]
            tag  = ''.join([key[part] for part in name.split('-') if part in key])
            print('retrieving {} as {}...'.format(name, tag))
            data = dataset.load(file).data
            S = struct(vars=get_vars_from(data), stats=get_stats_from(data), real=count_subjects_from(data), virtual=get_subject_number_from(data))
            if tag.startswith('go'):
                x, f, r = get_curves_from(data)
                S.set(x=x, f=f, r=r)
            auc  = AUC_of(data, 'before', 'bad', 1) if tag.startswith('go') else AUC_of(data)
            pst = show_(S.stats)
            S.set(auc=auc, pst=pst)
            if 'bio' in S.stats:
                auc = AUC_of(data, 'bio')
                S.set(bio=auc)
            R.set(**{tag:S})
            print()
    all_real, all_virtual = [], 0
    for _, data in R.tokens:
        all_real += data.real
        all_virtual += data.virtual
    print('{} virtual patients from {} real subjects assessed'.format(all_virtual, len(set(all_real))))
    return R

no_color = (0,0,0,0)
red =   (1,0,0)
green = colors.cnames['green']
hgreen = (0,.8,0)
blue =  (0,0,1)
grey =  (.3,.3,.3)
orange = colors.cnames['orange']
black = (0,0,0)
def show_curve_from(results, size=(5,3), lw=3, fsize=10):
    _ = lab.figure(figsize=size)
    ax = lab.subplot2grid((1,1), (0,0), 1, 1, frameon=False)
    x, l, d = results.gol.x, results.gol.r, results.god.r
    lab.plot(x, color=black, linewidth=lw)
    lab.plot(d, color=green, linewidth=lw)
    lab.plot(l, color=hgreen, linewidth=lw)
    lab.tick_params(
        axis =      'x',
        #rotation =  45,
        labelbottom = True
    )
    lab.xticks(
        [0,     31-2,   61-2,   92-2, 122-2,  153-2],
        ['March','April','May','June','July','August'],
        fontsize=fsize, weight='bold', ha='left')
    lab.show()

def cf_panel(A, B, fig=None, grid=(1,1), size=(5,5), pos=(0,0), cols=1, rows=1, color=red, lw=5, labels=('EHR','LR'), p=0, ref_p=.05, fsize=10):
    show = False
    if not fig: fig = lab.figure(figsize=size); show=True
    ax = lab.subplot2grid(grid, pos, colspan=cols, rowspan=rows, frameon=False)
    bp=lab.boxplot((A,B), widths=[.6,.6], positions=[0,1], showfliers=False, patch_artist=True, labels=labels)
    alphas = [.5 if np.median(A)<np.median(B) else .75, .75 if np.median(A)<np.median(B) else .5]
    if p>=ref_p: alphas = [.5]*2
    for n,item in enumerate(bp['boxes']): item.set(facecolor=colors.to_rgba(color,alphas[n]), edgecolor=no_color)
    for n,item in enumerate(bp['caps']): item.set(color=no_color)
    for n,item in enumerate(bp['whiskers']): item.set(color=colors.to_rgba(color,alphas[int(round(n/3))]), linewidth=lw)
    for item in bp['medians']: item.set(color=color, linewidth=lw)
    lab.tick_params(
        axis =  'x',
        which = 'both',
        bottom = False,
        top =    False,
        labelbottom = True
    )
    lab.xticks(fontsize=fsize, weight='bold')
    if show: lab.show()

def days_from(results, size=(10,10), yfsz=7):
    fig, grid = lab.figure(figsize=size), (1,3)
    colors = (green, orange, red)
    A, B = ['d','dc','dc'], ['l','lc','lc']
    for n, panel in enumerate(('go', 'mv', 'ic')):
        cf_panel(results.get(panel+A[n]).stats['before'],results.get(panel+B[n]).stats['before'], fig, grid, None, (0,n), 1,1, colors[n], 5, p=1)
        lab.yticks(fontsize=yfsz)
    lab.show()