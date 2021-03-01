## APPLICATION OF COP ON COVID19 DATA
# pylint: disable=no-self-argument, no-member

from core import some, no, value, struct
from cop import load_data, translation, patient, model
import datetime as dt

from matplotlib import pylab as lab
from scipy.signal import savgol_filter as smooth
from scipy import stats
import numpy as np

_MODE = None #select type of analysis when running module, set to None when importing!
#_MODE = 'proc-early_discharge'
#_MODE = 'proc-CoV2-early_discharge'

#_MODE = 'proc-CoV2-O2-support'
#_MODE = 'proc-CoV2-ICU'

proc_thresh = .5 # threshold for bad outcome selection, see fig. 1
sep = ';'
results_dir = "...somefolder"
retry_model = 3 # number of iterations

def make_model(by_vars, data, groups = {(0,):0, range(1,5):1}):
    ''' return initialized model ready to be trained through "calculate" '''
    T = translation(by_vars)
    pclass = patient(by_vars)
    pclass.load(data, T)
    return model(pclass, groups=groups)

DATE = {'year':(0,4), 'month':(4,6), 'day':(6,8)}
var_dir = '...somefolder/COVID19/raw/'
out_dir = '...somefolder/COVID19/data/'
REF, dcrC, labC, dCp, lCp = None, None, None, None, None

class subject:
    ''' service class to calculate curve '''
    admission, dismission = None, None
    def add(this, date):
        if this.admission == None: this.admission = date
        else:
            ref = this.admission
            days = (ref-date).days
            if days>0:
                this.dismission = this.admission
                this.admission = date
            elif days<0:
                if this.dismission == None: this.dismission = date
                else:
                    ref = this.dismission
                    days = (ref-date).days
                    if days<0: this.dismission = date

class curve:
    ''' service class to calculate curves '''
    def __init__(_, first, last):
        days = (last-first).days
        _.day = [[] for day in range(days)]
        _.days = days
        _.start = first
    def add(this, sub, prediction=0):
        at, span = (sub.admission-this.start).days, (sub.dismission-this.start).days
        span -= int(at + prediction)
        if span>0:
            while span:
                if at<len(this.day): this.day[at].append(sub)
                span -= 1
    def calculate(curve): return [len(set(day)) for day in curve.day]

def calculate_curve(results, calculate_curve_using_levels_up_to = 2, data_file='data.csv', check_file='raw_data.csv', tries=1000):
    ''' function to calculate occupation curve corrected by a model '''
    idc, datec = 0, 3
    subs, pcr = {}, None
    first, last = dt.date.today(), dt.date(2020,1,1)
    
    print('processing subjects', end='...')
    with open(out_dir+check_file, 'r') as file:
        found, lines = [], file.readlines()
        varz = [var.strip() for var in lines[0].split(sep)]
        outc = varz.index('XXX') # variable hidden for maintaining double blind
        for line in lines[1:]:
            fields = line.split(sep)
            ID = fields[idc].strip()
            if calculate_curve_using_levels_up_to>=int(fields[outc]):
                found.append(int(ID if ID != '' else 0))
        pcr = list(set(found))
    print('done')

    print('calculating stays', end='...')
    with open(var_dir+data_file, 'r') as file:
        for line in file.readlines()[1:]:
            fields = line.split(sep)
            ID = fields[idc].strip()
            if ID[0].isdigit():
                subid, dateraw = int(ID), fields[datec].strip()
                if subid in pcr:
                    if not subid in subs: subs[subid] = subject()
                    args = {dmy:int(dateraw[DATE[dmy][0]:DATE[dmy][1]]) for dmy in DATE}
                    sub, date = subs[subid], dt.date(**args)
                    if (date-first).days<0: first = date
                    if (date-last).days>0: last = date
                    sub.add(date)
    
    for sub in subs:
        if subs[sub].dismission == None:
            ref = subs[sub].admission
            try:
                subs[sub].dismission = dt.date(ref.year, ref.month, ref.day+1)
            except:
                subs[sub].dismission = dt.date(ref.year, ref.month+1, 1)
    print('done')

    print('creating reference curve', end='...')
    ref = curve(first, last)
    for sub in subs: ref.add(subs[sub])
    print('done')

    import scipy.stats as stats
    import numpy as np
    print('creating adjusted curves', end='...')
    normal = stats.shapiro(results.stats['before'])[1]
    days = np.average(results.stats['before']) if normal>0.05 else np.median(results.stats['before'])
    drng = days + (np.std(results.stats['before']) if normal>0.05 else stats.iqr(results.stats['before']))
    threshold = len(results.dismissed)/(len(results.dismissed)+len(results.in_room))
    def off(opt):
        if opt: return np.random.randint(0,drng)
        return days
    fixed_curves, ranged_curves = [],[]
    for cset,opt in ((fixed_curves,0),(ranged_curves,1)):
        test = tries
        while test:
            C = curve(first, last)
            for s in subs:
                check = np.random.rand()
                C.add(subs[s], off(opt) if check<=threshold else 0)
            cset.append(C.calculate())
            test -= 1
    print('done')

    return first, last, fixed_curves, ranged_curves, ref.calculate()

def calculate_models(retries, hub, calc_args, curve_args, message):
    ''' calculate models and curves '''
    results = []
    while retries:
        print('{} tries before saving...'.format(retries))
        print(message.format(M.DB.ids, M.DB.long_ids))
        hub.calculate(**calc_args)
        hub.best = hub.unstacked
        ir, dc = set([int(s.hn) for s in hub.in_room]), set([int(s.hn) for s in hub.dismissed])
        result = struct(vars=hub.vars_of, stats=hub.stats, in_room=ir, discharged=dc)
        first, last, fcs, rcs, rc = calculate_curve(hub, **curve_args)
        result.set(first=first, last=last, fixed_curves=fcs, ranged_curves=rcs, reference_curve=rc)
        hub.best, hub.unstacked = [], []
        results.append(result)
        retries -= 1
    return results

def calculate_models_only(retries, hub, calc_args, message):
    ''' calculate models '''
    results = []
    while retries:
        print('{} tries before saving...'.format(retries))
        print(message.format(M.DB.ids, M.DB.long_ids))
        hub.calculate(**calc_args)
        hub.best = hub.unstacked
        assessed = set([int(s.hn) for s in hub.assessed])
        result = struct(vars=hub.vars_of, stats=hub.stats, assessed=assessed)
        hub.best, hub.unstacked = [], []
        results.append(result)
        retries -= 1
    return results    

def get_curves_from(results):
    ''' function to average curves '''
    ref = results[0].reference_curve
    fix, rng = [], []
    for data in results:
        fix += data.fixed_curves
        rng += data.ranged_curves
    return ref, np.average(fix,0), np.average(rng,0)

def get_subjects_from(results, force_room_home=False):
    ''' get subject IDs from results '''
    room, home = [], []
    for data in results:
        room += list(data.in_room) if 'in_room' in data.sets else list(data.assessed)
        home += list(data.discharged) if 'discharged' in data.sets else []
    if not force_room_home and home == []: return room
    return room, home

def get_subject_number_from(results):
    ''' calculate assessed subjects' number '''
    room, home = get_subjects_from(results, True)
    return len(set(room+home))

def get_vars_from(results):
    ''' retrieves all variable weights from a series of results '''
    vars = {l:[] for l,_ in results[0].vars}
    for data in results:
        for l,_set in data.vars:
            vars[l] += _set
    return vars

def get_stats_from(results):
    ''' retrieves statistics from a series of results '''
    stat = {s:[] for s,_set in results[0].stats.items()}
    for data in results:
        for s,_set in data.stats.items():
            stat[s] += _set
    return stat

def show_(vars, force_avg=False):
    ''' show variable weights from results, calculates median or mean depending on distribution, user can force average '''
    results = struct()
    for label,_set in vars.items():
        normal = stats.shapiro(_set)[1] if len(_set)<5000 else stats.normaltest(_set)[1]
        value = np.average(_set) if normal>=0.05 or force_avg else np.median(_set)
        rng = np.std(_set) if normal>0.05 or force_avg else stats.iqr(_set)
        print('{}:\t{:.2f}±{:.3f} ({:.3f})'.format(label, value, rng, normal))
        results.set(**{label:struct(value=value, range=rng, normal=normal)})
    return results

def AUC_of(stats, target='before', comparator='good', _from=0):
    ''' calculates ROC-AUC of stats using a target and a comparator '''
    values = [v for source in stats for v in source.stats[target]]
    sens = 1-sum([int(v<_from) for v in values])/len(values)
    spec = np.average([s for source in stats for s in source.stats[comparator]])
    AUC = (sens+spec)/2
    print('sensitivity: {:.1%}'.format(sens))
    print('specificity: {:.1%}'.format(spec))
    print('ROC-AUC:     {:.1%}'.format(AUC))
    return struct(sens=sens, spec=spec, AUC=AUC)

def plot(*lines, smooth_days=0, degree=0):
    ''' plots average curve, possibility to smooth NOT used in figure 2 '''
    if smooth_days>0:
        if smooth_days>2:
            redraw = lines
            lines = []
            for line in redraw: lines.append(smooth(line, smooth_days, degree))
        lines = [np.array(line)/max(lines[0]) for line in lines]
    for line in lines:lab.plot(line)
    lab.show()

# USE MODE TO SELECT THE ANALYSIS
# err is the tolerated deviation from outcome balance at training, ideal error is 0.0
if _MODE == 'proc-early_discharge':
    M, err = make_model('...somefolder/vars.txt', 'somefolder/data.csv'), .25
    calcargs = dict(on=.7, virtual=int(M.DB.long_ids*.2), outcome_from=0, err=err, from_end=True, threshold=proc_thresh, skip=1, targets=[])
    curvargs = dict(calculate_curve_using_levels_up_to = 2, data_file='data.csv', check_file='raw_data.csv', tries=1000)
    results = calculate_models(retry_model, M, calcargs, curvargs, 'dcr-anticipate deterioration or recovery: {} subjects, {} in room')
    #value(results).save(results_dir+name)
elif _MODE == 'proc-CoV2-early_discharge':
    M, err = make_model('...somefolder/vars.txt', '...somefolder/data.csv'), .4
    calcargs = dict(on=.7, virtual=int(M.DB.long_ids/10), outcome_from=0, err=err, from_end=True, threshold=proc_thresh, skip=1, targets=[])
    curvargs = dict(calculate_curve_using_levels_up_to = 2, data_file='data.csv', check_file='raw_data.csv', tries=1000)
    results = calculate_models(retry_model, M, calcargs, curvargs, 'dcr-PCR+-anticipate deterioration or recovery: {} subjects, {} in room')
    #value(results).save(results_dir+name)
elif _MODE == 'proc-CoV2-O2-support':
    groups = {(0,1):0, range(2,6):1} #overrides original group cut, see fig. 1 for all possibilites
    M, err = make_model('...somefolder/vars.txt', '...somefolder/data.csv', groups), .25
    calcargs = dict(on=.6, virtual=int(M.DB.long_ids*.3), outcome_from=3, err=err, skip=0)
    results = calculate_models_only(retry_model, M, calcargs, 'dcr-predict OSi: {} subjects, {} in room')
    #value(results).save(results_dir+name)
elif _MODE == 'proc-CoV2-ICU':
    groups = {range(3):0, range(3,6):1} #overrides original group cut, see fig. 1 for all possibilites
    M, err = make_model('...somefolder/vars.txt', '...somefolder/dcr-lab-CoV2.csv', groups), .2
    calcargs = dict(on=.7, virtual=int(M.DB.long_ids/5), outcome_from=4, err=err, skip=0)
    results = calculate_models_only(retry_model, M, calcargs, 'dcr-predict ICU: {} subjects, {} in room')
    #value(results).save(results_dir+name)