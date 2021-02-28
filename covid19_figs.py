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

def load(key, ftype='.res'):
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

ehr_x = ['hn', 'start', 'bio', 'outcome']
ehr_j = ['PCR', 'PCR2']
lab_x = ['hn', 'outcome']
def table_of(data, vars, exclude=[], join=[], sep='\t'):
    print('raw table...')
    var = show_(data)
    print('done.\n\ncreating final data table...')
    _vars = {}
    with open(vars, 'r', encoding='latin1') as file:
        for line in file.readlines():
            fd, label = line.split(':')
            _, field = fd.split('\t')
            label = label.strip()
            field = field.strip()
            if field not in exclude and field in var.sets:
                _vars[field] = (label, var.get(field).value, var.get(field).range, var.get(field).normal>=0.05)
    if len(join)>1:
        base = _vars[join[0]]
        for field in join[1:]:
            _, val, rng, nrl = _vars[field]
            base = (base[0], base[1]+val, base[2]+rng, base[3] == nrl and nrl)
            _vars.pop(field)
        _vars[join[0]] = base
    all = sorted(_vars.values(), key=lambda x:x[1], reverse=True)
    lines = [sep.join(['variable','value','normal'])]
    lines += [sep.join([l[0], '{:.1%}±{:.1%}'.format(l[1],l[2]), 'yes' if l[3] else 'no']) for l in all]
    for line in lines: print(line)
    return lines

def compare(tag='comp', ftype='.res', exclude=['Tlim'], compare=1000, times=10, thresh=.05):
    R = struct()
    for file in listdir():
        if file.endswith(ftype) and tag in file:
            name = file[:-len(ftype)]
            data = dataset.load(file).data
            sets = [_set for _set in data.dcr.sets if _set not in exclude]
            for _set in sets:
                A, B = data.dcr.get(_set), data.lab.get(_set)
                repeat = times
                As, Bs, p = [], [], []
                while repeat:
                    a = [A[i] for i in np.random.randint(0,len(A), compare)]
                    b = [B[i] for i in np.random.randint(0,len(B), compare)]
                    an, bn = stats.shapiro(a)[1], stats.shapiro(b)[1]
                    normal = an>=thresh and bn>=thresh
                    p.append(stats.ttest_ind(a,b)[1] if normal else stats.mannwhitneyu(a,b)[1])
                    As, Bs = As+a, Bs+b
                    repeat -= 1
                an, bn = stats.normaltest(As)[1], stats.normaltest(Bs)[1]
                normal = an>=thresh and bn>=thresh
                Am = np.average(As) if normal else np.median(As)
                Bm = np.average(Bs) if normal else np.median(Bs)
                Ad = np.std(As) if normal else stats.iqr(As)
                Bd = np.std(Bs) if normal else stats.iqr(Bs)
                pr = sum([int(i<thresh) for i in p])/len(p)
                pm = np.average(p) if stats.shapiro(p)[1]>thresh else np.median(p)
                result = struct(Am=Am, Bm=Bm, Ad=Ad, Bd=Bd, pm=pm, pr=pr, p=p, As=As, Bs=Bs)
                print('{} (dcr vs. lab): {:.2f}±{:.2f} vs. {:.2f}±{:.2f}, p={:.3f} ({:.0%})'.format(
                    _set,               Am, Ad,             Bm, Bd,        pm,       pr))
                R.set(**{_set:result})
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
    fig = lab.figure(figsize=size)
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

def comp_fig(results, size=(15,7),
            labels={'alt':'ALT', 'ast':'AST', 'bnp':'BNP', 'crea': 'Cr.', 'crp':'CRP', 'dd':'D-dimer', 'fer':'Fer.', 'il6':'IL6', 'leu':'WBC', 'plat':'Plts', 'proca':'Procal.'},
            colors = (red, green, grey)):
    fig, grid = lab.figure(figsize=size), (8,33)
    for n,panel in enumerate(('crp', 'dd', 'il6')):
        cf_panel(results.get(panel).As, results.get(panel).Bs, fig, grid, None, (0,n*11), 10, 5, colors[n], 5, (labels[panel], '{:.0%}'.format(results.get(panel).pr)), results.get(panel).pm, fsize=15)
    for n,panel in enumerate(('alt','ast','bnp','crea','fer','leu','plat','proca')):
        cf_panel(results.get(panel).As, results.get(panel).Bs, fig, grid, None, (6,n*4), 3, 2, blue, 5, (labels[panel], '{:.0%}'.format(results.get(panel).pr)), results.get(panel).pm)
    lab.show()

key_loc = "G:/My Drive/COVID19/data/"
def correct_results_to(loc='translated/', key='key.csv', key_location=key_loc, ftype='.res', exclude=['dcr-lab-comp.res'], sep=';'):
    trans = {}
    with open(key_location+key, 'r') as file:
        for line in file.readlines():
            ori, translated = line.split(sep)
            trans[int(ori)] = int(translated)
    for file in listdir():
        if file not in exclude and file.endswith(ftype):
            res = dataset.load(file).data
            for data in res:
                fields = ['assessed'] if 'assessed' in data.sets else ['discharged', 'in_room']
                for field in fields:
                    new = [trans[sub] for sub in data.get(field)]
                    data.set(**{field:set(new)})
            dataset(res).save(loc+file)