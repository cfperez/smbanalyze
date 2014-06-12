from operator import itemgetter, attrgetter, methodcaller
from matplotlib.cbook import flatten
import cPickle as pickle

class numlist(list):
  '''
  Ordered (and sorted) container for experiments
  '''
  def __init__(self, iterable=[]):
    super(numlist, self).__init__(iterable)
    try:
      self.sort()
    except (AttributeError, KeyError, ValueError):
      pass

  def __getslice__(self, i, j):
    return self.__getitem__(slice(i, j))

  def __getitems__(self, key_tuple):
    return numlist(super(numlist, self).__getitem__(key) for key in key_tuple)

  def __getitem__(self, key):
    if isinstance(key, tuple):
      return self.__getitems__(key)
    out = super(numlist, self).__getitem__(key)
    if isinstance(out, list):
      return numlist(out)
    else:
      return out

  def next(self):
    try:
      return self.it.next()
    except (StopIteration, AttributeError):
      self.it = iter(self)
      return self.next()

  def __add__(self, other):
    return numlist(super(numlist, self).__add__(other))

  def __repr__(self):
    fmt = '\n '.join('{0}: {1}'.format(i, repr(item))
                      for i, item in enumerate(self))
    return '[' + fmt + ']'

  @classmethod
  def map(cls, func, arr, *more):
    '''Returns a numlist() using built-in map()
    '''
    return cls(map(func, arr, *more))

  def flatten(self):
    return numlist(flatten(self))

  def __getattr__(self, attr):
    return self.getattrs(attr)

  def getattrs(self, attr):
    return self.map(attrgetter(attr), self.has_attr(attr))

  def getitems(self, key):
    return self.map(itemgetter(key), self)

  def filter(self, condition):
    "Returns a numlist with experiments matching condition"
    cls = type(self)
    return cls(filter(condition, self))

  def get(self, name):
    "Get all attributes of experiments in numlist by name: mol.get('f')"
    cls = type(self)
    try:
      return cls.map(attrgetter(name), self)
    except AttributeError:
      return cls.map(itemgetter(name), self)
    else:
      raise AttributeError('Missing attribute {0} in a numlist element'.format(name))

  def set(self, name, value):
    for p in self:
      p[name] = value
    return self

  def call(self, action, *args, **kwargs):
    "Call function <action> with *args and **kwargs on all experiments"
    try:
      return map( methodcaller(action, *args, **kwargs), self )
    except AttributeError:
      raise AttributeError('Missing method {0} in a numlist element'.format(action))

  def save(self, filename):
    with open(filename, 'wb') as fh:
      pickle.dump(self, fh, protocol=2)
