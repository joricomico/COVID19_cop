# CORE version 0.1
## introducing ghost system (gcomparable, gprintable), which are either interfaces or
## abstract classes; since this is for reminding programmers, not to force them,
## the Abstract class is not employed.
## THIS is a working version that will be replaced soon with a more "functional" one...
## IMPORTANT NOTE: throughout the code, pythonic way is not employed,
## the most notable change is in class definitions, where methods definitions tend to be readable,
## e.g., in class struct, get(object, field) automatically explains itself, self is replaced by object.
## for this reason some pylint options are disabled.
# pylint: disable=no-self-argument, no-member

_DEBUG = 0

from inspect import isfunction, ismethod, isgeneratorfunction, isgenerator, isroutine
from inspect import isabstract, isclass, ismodule, istraceback, isframe, iscode, isbuiltin
from inspect import ismethoddescriptor, isdatadescriptor, isgetsetdescriptor, ismemberdescriptor
from inspect import isawaitable, iscoroutinefunction, iscoroutine

from datetime import timedelta as _time
from datetime import datetime
from collections.abc import Iterable as iterable

from pickle import dump, load

def some(field): 
    ''' returns True if value is not None or pointing to an empty set; therefore 0, True and False return True '''
    return field is not None and field != [] and field != {} and field != () and field != ''
def no(field):
    ''' returns False if value is not None or pointing to an empty set; therefore 0, True and False return False '''
    return not some(field)

class clonable:
    ''' inherit from this to implement object cloning and metadata '''
    def __init__(clonable, **sets): clonable.__dict__.update(sets)
    def _clonable(get): return get.__dict__.copy()
    def _meta(data): return data._clonable()
    def clone(_): return type(_)(**_._clonable())
    def set(object, **fields):
        for field in fields: setattr(object, field, fields[field])
    @property
    def sets(of): return sorted(list(set(dir(of)) - set(dir(type(of)))))

class gcomparable:
    ''' ghost: inherit from this to implement object comparison by internal dictionary '''
    def _compare(a, b): 
        if type(a) != type(b): return False
        if a.__dict__ == b.__dict__: return True
        return False
    def __eq__(a, b): return a._compare(b)

class gprintable:
    ''' ghost: inherit from this to port some terminal printing functions '''
    _lines_ = 31
    _chars_ = 13
    _ellipsis_ = '...'
    def _repr(my, value):
        _type = ''.join(''.join(str(type(value)).split('class ')).split("'"))
        _value = '{}'.format(value)
        if len(_value)>my._chars_:
            show = int(my._chars_/2)
            _value = _value[:show]+my._ellipsis_+_value[-show:]
        return '{} {}'.format(_type, _value)

class struct(clonable, gcomparable, gprintable):
    ''' main objecto to store data '''
    @staticmethod
    def _from(data):
        if value(data).inherits(struct): return struct(**data._clonable())
        elif hasattr(data, '__dict__'): return struct(**data.__dict__)
        return value(data)
    def _default(field, name, value):
        try: return getattr(field, name)
        except: setattr(field, name, value)
        return value
    def all(object, *fields): return [getattr(object, field) for field in fields if field in object.__dict__]
    def get(object, field):
        if field in object.sets: return getattr(object, field)
        return None
    def _check(these, args, by=lambda x:x):
        def match(this, item, value): return item in this.sets and by(this.get(item)) == value
        return all([match(these, _, args[_]) for _ in args])
    def clear(_, *fields):
        if no(fields): fields = _.sets
        for field in [field for field in fields if hasattr(_,field) and not ismethod(getattr(_, field))]: delattr(_, field)
    @property
    def tokens(_): return ((k,_.get(k)) for k in _.sets)
    def __repr__(self):
        if not hasattr(self, '_preprint'): return struct(_preprint='', _lines=self._lines_, data=self).__repr__()
        pre, repr = self._preprint, ''
        for n,i in ni(self.data):
            if self._lines == 0: break
            else: self._lines -= 1
            repr += pre+'{}: '.format(n)
            if issubclass(type(i), struct): repr += '\n'+struct(_preprint=pre+'\t', _lines=self._lines, data = i).__repr__()
            else: repr += self._repr(i)
            repr += '\n'
        return repr

class recordable(clonable):
    ''' inherit from this to enable load, save and other methods for storage '''
    @staticmethod
    def load(filename):
        with open(filename, 'rb') as file: return load(file)
    def _to_base(_, value): return value
    def _parse_to_base(_):
        clonable = _._clonable()
        for field in clonable:
            if issubclass(type(field), recordable): clonable[field] = clonable[field]._parse_to_base()
            else: clonable[field] = _._to_base(clonable[field])
        return type(_)(**clonable)
    def _predump(_): pass
    def save(data, filename, to_base=False):
        if to_base: _ = data._parse_to_base()
        data._predump()
        with open(filename, 'wb') as file: dump(data, file)
        
class value(recordable, gcomparable, gprintable):
    ''' main class for reflection and storage '''
    data = None
    _check = dict(
            isfunction=isfunction, ismethod=ismethod, isgeneratorfunction=isgeneratorfunction, isgenerator=isgenerator, isroutine=isroutine,
            isabstract=isabstract, isclass=isclass, ismodule=ismodule, istraceback=istraceback, isframe=isframe, iscode=iscode, isbuiltin=isbuiltin,
            ismethoddescriptor=ismethoddescriptor, isdatadescriptor=isdatadescriptor, isgetsetdescriptor=isgetsetdescriptor, ismemberdescriptor=ismemberdescriptor,
            isawaitable=isawaitable, iscoroutinefunction=iscoroutinefunction, iscoroutine=iscoroutine
                   )
    def __init__(this, token, **meta):
        this.data = token
        this.__dict__.update({k:v(token) for k,v in this._check.items()})
        super().__init__(**meta)
    @property
    def type(_): return type(_.data)
    def inherits(_, *types): return issubclass(_.type, types)
    @property
    def isstruct(_): return _.inherits(struct)
    @property
    def isbaseiterable(_): return _.inherits(tuple, list, dict, set) or _.isgenerator or _.isgeneratorfunction
    @property
    def isiterable(_): return isinstance(_.data, iterable) and _.type is not str
    def _clone_iterable(_):
        if _.inherits(dict): return _.data.copy()
        elif _.isgenerator or _.isgeneratorfunction: return (i for i in list(_.data))
        else: return type(_.data)(list(_.data)[:])
    def _clonable(_): return {k:v for k,v in _.__dict__.items() if k not in _._check}
    def _meta(data): return {k:v for k,v in data._clonable().items() if k != 'data'}
    def clone(_):
        data = _.data
        if _.isiterable: data = _._clone_iterable()
        elif _.inherits(clonable): data = _.data.clone()
        return type(_)(data)
    def __enter__(self): self._instance = self; return self
    def __exit__(self, type, value, traceback): self._instance = None
    def __repr__(self):
        if not hasattr(self, '_preprint'): return value(self.data, _preprint='', _lines=value(value._lines_)).__repr__()
        if self.isbaseiterable:
            pre, repr = self._preprint, ''
            for n,i in ni(self.data):
                if self._lines.data == 0: break
                else: self._lines.data -= 1
                index, item = str(n), i
                if self.inherits(dict): index += ' ({})'.format(str(i)); item = self.data[i]
                repr += pre+'{}: '.format(index)
                next = value(item, _preprint=pre+'\t', _lines=self._lines)
                if next.isiterable: repr += '\n'
                repr += next.__repr__()
                repr += '\n'
            return repr
        elif self.inherits(clonable): return value(self.data._clonable(), _preprint=self._preprint, _lines=self._lines).__repr__()
        else: return self._repr(self.data)
this = value

def meta(data):
    if this(data).inherits(clonable): return data._meta()
    return struct._from(data)._meta()
def get(opt, key, default=None, share=False):
    if key in opt:
        if not share: return opt.pop(key)
        return opt[key]
    return default

def ni(list):
    if this(list).isiterable:
        for n,i in enumerate(list): yield n,i
    elif this(list).isstruct:
        for n,i in list.tokens: yield n,i
    else: yield None, list

class at(struct):
    DAY, HOUR, MIN = 86400, 3600, 60
    def __init__(_, dtime=None, **sets):
        super().__init__(**sets)
        if some(dtime) and issubclass(type(dtime), _time): _._time = dtime
        else:
            d,h,m,s,ms = _._default('d',0), _._default('h',0), _._default('m',0), _._default('s',0), _._default('ms',0)
            if not any([d,h,m,s,ms]): now=datetime.now(); _._time = now-datetime(now.year, now.month, now.day)
            else: _._time = _time(days=d, hours=h, minutes=m, seconds=s, milliseconds=ms)
        _.clear('d','h','m','s','ms')
    def __sub__(_, dtime):
        of=type(dtime); sets=_._clonable()
        if issubclass(of, _time): return at(_._time-dtime, **sets)
        elif issubclass(of, at): sets.update(dtime._clonable()); return at(_._time-dtime._time, **sets)
    def __add__(_, dtime):
        of=type(dtime); sets=_._clonable()
        if issubclass(of, _time): return at(_._time+dtime, **sets)
        elif issubclass(of, at): sets.update(dtime._clonable()); return at(_._time+dtime._time, **sets)
    def __str__(_): return str(_._time)
    @property
    def seconds(_): return _._time.seconds
    @property
    def S(_): return _.seconds
    @property
    def minutes(_): return _._time.seconds/60
    @property
    def M(_): return _.minutes
    @property
    def hours(_): return _.minutes/60
    @property
    def H(_): return _.hours
    @property
    def days(_): return _._time.days
    @property
    def D(_): return _.days
    @staticmethod
    def zero(): return at(_time())

if _DEBUG:
    print(some(None), some(0), some(True), some(False), some([]), some(tuple()))
    print(no(None), no(0), no(True), no(False), no([]), no(tuple()))