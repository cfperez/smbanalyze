import collections

from numpy import all, where, asarray, sign, ndarray

from fileIO import load, toSettings


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
  def fromFile(cls, filename):
    meta, data = load(filename, comments=toSettings)
    meta[cls.__name__+'_filename'] = filename
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
      raise AttributeError("{0} has no attribute '{1}'".format(
          self.__class__.__name__, attr))

  @property
  def T(self):
    return self.data.T

  def __iter__(self):
    return iter(self.T)

  def __getitem__(self, key):
    return self.data[key].T

  def __repr__(self):
    return repr(self.data)

class TrapData(AbstractData):
  _fields = TrapData_fields

  def maskFromLimits(self, x, f, limits=(0,-1)):
    start, stop = limits
    ext_fit,f_fit = self.ext[start:stop], self.f[start:stop]

    if f is None: f=[max(f_fit)]
    try:
      min_f, max_f = f
    except TypeError:
      min_f, max_f = min(f_fit), f
    except ValueError:
      min_f, max_f = min(f_fit), f[0]
    if max_f>max(f_fit): max_f=max(f_fit)

    if x is None: x=[min(ext_fit)]
    try:
      min_ext, max_ext = x
    except TypeError:
      min_ext, max_ext = x, ext_fit[min(where(f_fit>=max_f)[0])-1]
    except ValueError:
      min_ext, max_ext = x[0], ext_fit[min(where(f_fit>=max_f)[0])-1]

    between = lambda s,a,b: (s>=a) & (s<=b)
    mask = between(ext_fit, min_ext, max_ext)

    return mask

class FretData(AbstractData):
  _fields = FretData_fields


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
