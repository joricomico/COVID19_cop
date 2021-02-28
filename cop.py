# COP version 0.0
# pylint: disable=no-self-argument, no-member
 
from core import struct, some, no, ni

from numpy import array, average, std
from numpy.random import rand, randint, seed

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, KFold, RepeatedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from sklearn.pipeline import Pipeline

import datetime as dt

_DEBUG = 0

def load_data(file, seps='\t,;'):
    n_seps = len(seps)
    data = []
    try:
        raw = open(file, 'r', encoding='utf8')
        lines = raw.readlines()
    except:
        raw = open(file, 'r')
        lines = raw.readlines()
    for line in lines:
        tokens = []
        for n in range(n_seps):
            test = line.split(seps[n])
            if len(test)>len(tokens): tokens = test
        if tokens[-1] == '\n': tokens.pop(-1)
        data.append(tokens)
    return data

class translation(struct):
    zero = 0.0
    _c2d = str.maketrans(',','.')
    def __init__(tclass, _from, at={'tags':1, 'types':2}):
        if some(_from) and some(at):
            data = load_data(_from)
            tags = data[at['tags']]
            types = data[at['types']]
            if len(tags) != len(types): raise ValueError('Tag Type mismatch:', len(types)-len(tags))
            super().__init__(tags=[t.strip() for t in tags], types=[t.strip() for t in types], items={})
    def translate(_, line, zero=None):
        if no(zero): zero = _.zero
        def _translate(this, value):
            value, tag = value.strip(), _.tags[this]
            if no(value): return zero
            _type, items = _.types[this], _.items
            if _type == 'n':
                try: return float(value.translate(_._c2d))
                except: return zero
            elif _type.startswith('d'):
                sep = _type[1]
                rawdate = value.split(sep)
                ref = dt.date.today()
                args = {}
                for n,in_type in enumerate(_type[2:]):
                    if in_type == 'd': args['day'] = int(rawdate[n])
                    elif in_type == 'm': args['month'] = int(rawdate[n])
                    else: args['year'] = int(rawdate[n])
                date = dt.date(**args)
                return float((ref-date).days)
            elif _type == 'l':
                if tag not in items: items[tag] = {}
                if value not in items[tag]: items[tag].update({value: len(items[tag])})
                return float(items[tag][value])
            elif _type.startswith('-'): return float(value.strip(_type[1:]))
            elif _type.endswith('-'):
                if value.strip() == _type.strip('- '): return 1.
                return zero
            else:
                if value.strip() == _type.strip(): return zero
                return 1.
        if len(line) == len(_.tags): return {tag: _translate(n,line[n]) for n,tag in ni(_.tags)}
        return {}

class patient(struct):
    ''' main structure to handle patient database; it's both root by statically using "states", and branch by using "prev" and "next" '''
    DEF = None
    states = []
    prev, next = None, None
    zero = 'prev'
    ID = None
    @staticmethod
    def reset(): patient.states = []
    def __init__(pclass, _from=None, line=1, args=None):
        if some(_from):
            args = load_data(_from)[line]
            patient.DEF = {arg.strip():0 for arg in args}
            patient._lines_ = len(patient.DEF)
        elif some(pclass.DEF):
            super().__init__(**pclass.DEF)
            pclass.states.append(pclass)
            pclass.set(**args)
    @property
    def len(_):
        if no(_.DEF): raise Warning('Definition not set')
        return len(_.DEF)
    def link(this, screening):
        this.next = screening
        screening.prev = this
    def _count_id(_, n=0):
        if some(_.ID): return _.ID, n
        return _.prev._count_id(n+1)
    @property
    def id(_): return _._count_id()
    def load(data, file, translator, id_from='hn', prev_at=0, _from=1):
        last, lines = None, load_data(file)[_from:]
        for line in lines:
            raw, id = translator.translate(line), id_from
            new = patient(args=raw)
            if raw[id] == prev_at and some(last):                
                reset = translator.translate(line, data.zero)
                for field in raw:
                    if reset[field] == data.zero: reset[field] = last.get(field)
                new.set(**reset)
                last.link(new)
            else: new.ID = raw[id]
            last = new
    def select(_, x_by='day', fx=lambda x:True, exclude=['hn'], y_by='outcome', groups={(0,1,2):0, (3,4,5):1}, prediction=0):
        def search_outcome_of(state, _from, this):
            check = lambda x,y:x>y
            if _from<0: check = lambda x,y:x<y
            search = state
            while search.prev:
                out = search.get(y_by)
                if check(out, this): this = out
                search = search.prev
            search = state
            while search.next:
                out = search.get(y_by)
                if check(out, this): this = out
                search = search.next
            return this
        X, Y, off = [], [], exclude + [y_by]
        for state in _.states:
            if fx(state.get(x_by)):
                out, y = state.get(y_by), None
                if prediction: out = search_outcome_of(state, prediction, out)
                for this in groups:
                    if out in this: y = groups[this]
                if no(y): continue
                x = [state.get(var) for var in state.sets if var not in off]
                X.append(x)
                Y.append(y)
        return X,Y
    @property
    def ids(_): return len([state for state in _.states if state.ID is not None])
    @property
    def long_ids(_): return len([state for state in _.states if state.ID and state.next])

def idsets_of(db):
    ''' function to yield single patients, patients are a series of states, whose only the first has valid ID and represent the subject '''
    _set = [db.states[0]]
    for state in db.states[1:]:
        if state.ID:
            yield _set
            _set = [state]
        else: _set.append(state)
    yield _set

class model(struct):
    ''' model structure: this will be updated to work online, so far it provides analysis and model construction '''
    k, err = 6, .15
    var_mul = 4
    ideal_y = .5
    classifier=RandomForestClassifier
    scaler=StandardScaler
    random_state = 31
    cargs = {}
    auto_exclude = ['hn', 'start']
    stack, unstacked = False, []
    _ref, from_end = 0, False
    train_on_head = False
    @staticmethod
    def score_of(x,y, margs={'n_estimators':100}):
        _model, scaler = model.classifier, model.scaler
        random_state = model.random_state
        if no(random_state): random_state = randint(0xFFFFFFFF)
        C = Pipeline([('scaler', scaler()), ('model', _model(**margs))])
        cv = KFold(model.k, shuffle=True, random_state=random_state)
        return cross_val_score(C, x, y, cv=cv)
    def __init__(M, database, groups={range(0,2):0, range(2,5):1}, by='outcome', prediction=1, exclude=[], **opts):
        super().__init__(**opts)
        exclude = list(set(exclude + M.auto_exclude))
        features = [feat for feat in database.states[0].sets if feat not in exclude and feat != by]
        M.set(_db=database, _selargs={'x_by':'day', 'groups':groups, 'y_by':by, 'exclude':exclude, 'prediction':prediction}, features=features, groups=groups, outcome=by, best=[])
    def train(_, on=.5, max_ids=0, err=.15, **opts):
        _.set(**opts)
        if _.random_state: seed(_.random_state)
        if max_ids == 0: max_ids = _._db.ids
        times, nest = _.k, _._db.len*_.var_mul
        DB, states = _._db, _._db.states
        all_states = states
        if _.train_on_head: states = [state for state in states if state.ID]
        sort, models = int(on*len(states)), []
        miny, maxy = _.ideal_y-err, _.ideal_y+err
        scores = []
        while times:
            sel = set(randint(0, len(states), sort))
            A,B = [state for n,state in ni(states) if n in sel], [state for n,state in ni(states) if n not in sel]
            DB.states = A
            if not _.train_on_head and DB.ids > max_ids: continue
            x,y = DB.select(**_._selargs)
            ay = average(y)
            if not _.train_on_head and ay<miny or ay>maxy: continue
            X,Y = array(x), array(y)
            S = _.scaler()
            Xs = S.fit_transform(X)
            M = _.classifier(n_estimators=nest, **_.cargs)
            M.fit(Xs,Y)
            DB.states = B
            x,y = DB.select(**_._selargs)
            Xt,Yt = array(x), array(y)
            Xts = S.transform(Xt)
            Yp = M.predict(Xts)
            score = accuracy_score(Yt, Yp)
            print('{:.3f}'.format(score), end='\t')
            models.append((M,S)); scores.append(score)
            times -= 1
        print('| {:.3f}'.format(average(scores)))
        ms = zip(models, scores)
        s_ms = sorted(ms, key=lambda x:x[1], reverse=True)
        DB.states = all_states
        _.best += s_ms
    def _calculate(_, on, by, targets, outcome_from, model, models, good, bad, threshold, _to=0):
        ''' partial method, retrieves the number of predicted subjects, this will be edited for online prediction '''
        start = _to
        DB, _._states = _._db, _._db.states
        train, test, p = [], [], [int(r>=on) for r in rand(_._db.ids)]
        for n, idset in ni(idsets_of(DB)):
            if p[n]: test.append(idset)
            else: train += idset
        DB.states = train
        _.train(on=model, err=_.err)
        if no(models): models = len(_.best)
        for idset in test:
            DB.states = idset
            (x,y), ym = DB.select(**_._selargs), []
            for ((M,S),s) in _.best[:models]:
                try:
                    X = array(x)
                    Xs = S.transform(X)
                    Yp = M.predict(Xs)
                    ym.append(Yp)
                except:
                    pass
            if len(ym) == 0: continue
            Y = array(ym).T
            p = average(Y)
            if not y[0]: p = 1-p
            if y[0]==good[0]: _.stats[good[1]].append(p)
            else: _.stats[bad[1]].append(p)
            if not _.from_end: _.assessed.append(idset[0])
            if len(idset)==1 or y[0]==_._ref or len(Y)!=len(idset): continue
            ref, target_map = None, {tag:None for tag in targets+['before']}
            #assert(len(Y)==len(idset))
            for n,state in ni(idset):
                yn = average(Y[n])
                if y[0]==0: yn = 1-yn
                if no(state.next) or yn>=threshold:
                    ref = state.get(by)
                    break
            if _.from_end:
                for state in idset:
                    if no(state.next):
                        anticipation = state.get(by)-ref
                        if anticipation:
                            target_map['before'] = anticipation
                            _.dismissed.append(idset[0])
                        else: _.in_room.append(idset[0])
            else:
                for state in idset:
                    if no(target_map['before']) and outcome_from<=state.get(_.outcome):
                        target_map['before'] = state.get(by)-ref
                        break
            for target in targets:
                for state in idset:
                    if no(target_map[target]) and state.get(target):
                        target_map[target] = state.get(by)-ref
                        break
            for tag in target_map:
                if target_map[tag] is not None: _.stats[tag].append(target_map[tag])
            _to+=1
        print('from {} to {}...'.format(start, _to))
        for tag in _.stats:
            print(tag+'\t{:.3f}'.format(average(_.stats[tag]) if len(_.stats[tag])>0 else 0))
        DB.states = _._states
        return _to
    #this is a batch method to run a series of iterations (_calculate)
    #as in other methods, pythonic style is not used, self is substituted by more human readable term:
    #"calculate stats on 50% (balance) by day (variable), targeting bio (additional target), ...
    # (considering) outcome_from 4, "
    def calculate(stats, on=.5, by='day', targets=['bio'], outcome_from=4, model=.9, virtual=None, models=None, good=(0,'good'), bad=(1,'bad'), skip=0, threshold=.75, **opts):
        stats.set(**opts)
        stats._ref = skip
        stats.set(stats = {tag:[] for tag in [good[1], bad[1], 'before']+targets})
        if stats.from_end: stats.set(in_room=[], dismissed=[])
        else: stats.set(assessed=[])
        n=stats._calculate(on, by, targets, outcome_from, model, models, good, bad, threshold)
        while some(virtual) and n<virtual:
            n=stats._calculate(on, by, targets, outcome_from, model, models, good, bad, threshold, n)
            if not stats.stack: stats.unstacked += stats.best; stats.best = []
    def reset(_):
        _._db.states = _._states
        _.best = []
    def sorted(_, model=0):
        (M,S), score = _.best[model]
        fi = zip(_.features, list(M.feature_importances_))
        s_fi = sorted(fi, key=lambda x:x[1], reverse=True)
        return s_fi, score
    def swap_stack(_):
        if len(_.best)==0: _.best = _.unstacked
        else: _.best = []
    @property
    def vars_of(model, start=0, stop=None):
        vari, s = model.sorted(start)
        index = {var:[i] for var,i in vari}
        if no(stop): stop = len(model.best)
        for n in range(start+1,stop):
            vari, s = model.sorted(n)
            for var,i in vari: index[var].append(i)
        s_vari = sorted(index.items(), key=lambda x:average(x[1]), reverse=True)
        return s_vari
    def show(_):
        for var,i in _.vars_of:
            print('{}\t{:.3f}'.format(var,average(i)))
    @property
    def DB(_): return _._db

if _DEBUG == 1:
    #HOW TO LOAD VARIABLES AND MODEL DATA
    T = translation('...somefolder/vars.txt')
    pclass = patient('...somefolder/vars.txt')
    pclass.load('...somefolder/data.csv', T)