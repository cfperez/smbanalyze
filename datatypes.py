import collections

from numpy import all, where, asarray, sign, ndarray, vstack
from operator import isSequenceType

from fileIO import load, toSettings, fileIOError


TYPE_CHECKING = 'STRICT'

FretData_fields = ('time','donor','acceptor','fret')
TrapData_fields = ('ext','f','sep')

class AbstractData(object):
  def __init__(self, data, **meta):
    self.data = data
    self.metadata = meta

  @classmethod
  def fromObject(cls, obj):
    try:
      return cls(obj.data, **obj.metadata)
    except AttributeError:
      raise ValueError('Constructor method only takes object with AbstractData interface')

  @classmethod
  def name(cls):
    return cls.__name__.lower()

  @classmethod
  def fromFile(cls, filename):
    try:
      meta, data = load(filename, comments=toSettings)
    except fileIOError as e:
      if not e.isError:
        print e.strerror
        data = load(filename)
        meta = {}
      else:
        raise
    else:
      meta[cls.name()+'_filename'] = filename
    finally:
      me = cls(data, **meta)
      return me

  @classmethod
  def fromFields(cls, *args, **meta):
    if len(args) != len(cls._fields):
      raise ValueError(
        "Number of arguments to {0} must match _fields".format(cls.__name__))
    return cls(asarray(args).T, **meta)

  def copy(self):
    return self.__class__.fromObject(self)

  def __len__(self):
    return len(self.data)

  def __eq__(self, other):
    return all(self.data==other.data) and self.metadata==other.metadata

  def __ne__(self, other):
    return not self==other

  def __add__(self, other):
    meta = self.metadata.copy()
    meta.update(other.metadata)
    return type(self)( vstack((self.data,other.data)) )


  def at(self, **kwargs):
    if len(kwargs)>1:
      raise ValueError('Use only one keyword representing a field')

    field, value = kwargs.popitem()
    if field not in self._fields:
      raise ValueError('Keyword argument must be a field')

    index = search_monotonic(getattr(self,field), value)
    return self[index]

  def __getattr__(self, attr):
    if attr in self._fields:
      return self.data[:,self._fields.index(attr)]
    else:
      try:
        return getattr(super(AbstractData, self), attr)
      except AttributeError:
        raise AttributeError("{0} has no attribute '{1}'".format(
            self.__class__.__name__, attr))

  @property
  def T(self):
    return self.data.T

  def __iter__(self):
    return iter(self.T)

  def __getitem__(self, key):
    return type(self)( self.data[key].view(), **self.metadata )

  def __repr__(self):
    return repr(self.data)

  @classmethod
  def _normalizeLimits(cls, limits, min_max, assume_max_limit=True):
    if limits is None:
      return min_max
    elif not isSequenceType(limits):
      limits = [limits]
    if len(limits) == 2:
      return limits
    elif len(limits) == 1:
      if assume_max_limit:
        return min_max[0], limits[0]
      else:
        return limits[0], min_max[1]

class TrapData(AbstractData):
  _fields = TrapData_fields

  def maskFromLimits(self, x, f, limits=(0,-1)):
    start, stop = limits
    ext_fit,f_fit = self.ext[start:stop], self.f[start:stop]

    min_f, max_f = TrapData._normalizeLimits( f,
                      min_max=(min(f_fit), max(f_fit)),
                      assume_max_limit=True
                    )

    min_ext, max_ext = TrapData._normalizeLimits( x,
                      min_max=(min(ext_fit), max(ext_fit)),
                      assume_max_limit=False
                    )

    between = lambda s,a,b: (s>=a) & (s<=b)
    return between(ext_fit, min_ext, max_ext) & between(f_fit, min_f, max_f)

  def select(self, x=None, f=None, limits=(0,-1)):
    return self[self.maskFromLimits(x, f, limits)]

  @property
  def fec(self):
    x,f,s = self.T
    return x,f

class FretData(AbstractData):
  _fields = FretData_fields

  def maskFromLimits(self, time, limits=(0,-1)):
    return 

  def select(self, time=None):
    return self[self.maskFromLimits(time)]


def search_monotonic(ar, value):
  shifted = ar-value
  if value <= ar[0]: return 0
  elif value >= ar[-1]: return -1

  start_sign = sign(shifted[0])
  for n, (current,last) in enumerate(and_prev(shifted)):
    if sign(current) != start_sign:
      return n if abs(current)<abs(last) else min(n-1, 0)
  return -1

def and_prev(iterable, default=None):
  last = default
  for x in iterable:
    yield x, last
    last = x

def _hasData(datatype):
  fields = getattr(datatype,'_fields',datatype)
  if TYPE_CHECKING=='STRICT':
    return lambda obj: all(map(hasattr,[obj]*len(fields),fields))
  else:
    return lambda obj: len(fields)==len(obj)

hasFretData = _hasData(FretData)
hasTrapData = _hasData(TrapData)
hasTrapFretData = lambda x: hasFretData(x) and hasTrapData(x)
