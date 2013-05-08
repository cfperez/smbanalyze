from functools import wraps
from inspect import getargspec
from itertools import izip_longest, takewhile

import numpy as np

def toNum(s):
  if s is None:
    return None
  try:
    return int(s)
  except ValueError:
    return float(s)

def toInt(s):
  s = toNum(s)
  return int(s) if s is not None else None

def isInt(s):
  try:
    return s == str(toInt(s))
  except ValueError:
    return False

def negate(it):
  return [-i for i in it]

def broadcast(f):
  @wraps(f)
  def _f(X,*args,**kwargs):
    if np.iterable(X):
      return np.array( map(lambda x: f(x,*args,**kwargs),X) )
    else:
      return f(X,*args,**kwargs)
  _f.arglist = getargspec(f).args
  _f.broadcast = True
  return _f

def fix_args(f, **fixed):
  ind_var = [fixed.pop('independent_var', 'x')]
  f_args = getattr(f, 'arglist', None) or getargspec(f).args
  pos_args = filter(lambda x: x not in fixed, f_args)

  if getattr(f, 'broadcast', False):
    @wraps(f)
    def _f(*args, **kwargs):
      newkwargs = kwargs.copy()
      newkwargs.update(zip(pos_args,args), **fixed)
      first = newkwargs.pop(pos_args[0])
      return f(first, **newkwargs)
  else:
    @wraps(f)
    def _f(*args, **kwargs):
      newkwargs = kwargs.copy()
      newkwargs.update(zip(pos_args,args), **fixed)
      return f(**newkwargs)

  return _f

def makeMatchStrFromArgs(*globs, **options):
  globs = list(globs)
  last = globs[-1]
  if isInt(last):
    globs[-1] = '_'+last
  if options.get('re_match', True):
    anychar = '.*'
    endmatch = r'(_|$)'
  else:
    anychar = '*'
    endmatch = ''
  return r'{any}{pattern}{end}'.format(pattern=anychar.join(globs), any=anychar, end=endmatch)

## {{{ http://code.activestate.com/recipes/576693/ (r9)
# Backport of OrderedDict() class that runs on Python 2.4, 2.5, 2.6, 2.7 and pypy.
# Passes Python2.7's test suite and incorporates all the latest updates.

try:
    from thread import get_ident as _get_ident
except ImportError:
    from dummy_thread import get_ident as _get_ident

try:
    from _abcoll import KeysView, ValuesView, ItemsView
except ImportError:
    pass


class OrderedDict(dict):
    'Dictionary that remembers insertion order'
    # An inherited dict maps keys to values.
    # The inherited dict provides __getitem__, __len__, __contains__, and get.
    # The remaining methods are order-aware.
    # Big-O running times for all methods are the same as for regular dictionaries.

    # The internal self.__map dictionary maps keys to links in a doubly linked list.
    # The circular doubly linked list starts and ends with a sentinel element.
    # The sentinel element never gets deleted (this simplifies the algorithm).
    # Each link is stored as a list of length three:  [PREV, NEXT, KEY].

    def __init__(self, *args, **kwds):
        '''Initialize an ordered dictionary.  Signature is the same as for
        regular dictionaries, but keyword arguments are not recommended
        because their insertion order is arbitrary.

        '''
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__root
        except AttributeError:
            self.__root = root = []                     # sentinel node
            root[:] = [root, root, None]
            self.__map = {}
        self.__update(*args, **kwds)

    def __setitem__(self, key, value, dict_setitem=dict.__setitem__):
        'od.__setitem__(i, y) <==> od[i]=y'
        # Setting a new item creates a new link which goes at the end of the linked
        # list, and the inherited dictionary is updated with the new key/value pair.
        if key not in self:
            root = self.__root
            last = root[0]
            last[1] = root[0] = self.__map[key] = [last, root, key]
        dict_setitem(self, key, value)

    def __delitem__(self, key, dict_delitem=dict.__delitem__):
        'od.__delitem__(y) <==> del od[y]'
        # Deleting an existing item uses self.__map to find the link which is
        # then removed by updating the links in the predecessor and successor nodes.
        dict_delitem(self, key)
        link_prev, link_next, key = self.__map.pop(key)
        link_prev[1] = link_next
        link_next[0] = link_prev

    def __iter__(self):
        'od.__iter__() <==> iter(od)'
        root = self.__root
        curr = root[1]
        while curr is not root:
            yield curr[2]
            curr = curr[1]

    def __reversed__(self):
        'od.__reversed__() <==> reversed(od)'
        root = self.__root
        curr = root[0]
        while curr is not root:
            yield curr[2]
            curr = curr[0]

    def clear(self):
        'od.clear() -> None.  Remove all items from od.'
        try:
            for node in self.__map.itervalues():
                del node[:]
            root = self.__root
            root[:] = [root, root, None]
            self.__map.clear()
        except AttributeError:
            pass
        dict.clear(self)

    def popitem(self, last=True):
        '''od.popitem() -> (k, v), return and remove a (key, value) pair.
        Pairs are returned in LIFO order if last is true or FIFO order if false.

        '''
        if not self:
            raise KeyError('dictionary is empty')
        root = self.__root
        if last:
            link = root[0]
            link_prev = link[0]
            link_prev[1] = root
            root[0] = link_prev
        else:
            link = root[1]
            link_next = link[1]
            root[1] = link_next
            link_next[0] = root
        key = link[2]
        del self.__map[key]
        value = dict.pop(self, key)
        return key, value

    # -- the following methods do not depend on the internal structure --

    def keys(self):
        'od.keys() -> list of keys in od'
        return list(self)

    def values(self):
        'od.values() -> list of values in od'
        return [self[key] for key in self]

    def items(self):
        'od.items() -> list of (key, value) pairs in od'
        return [(key, self[key]) for key in self]

    def iterkeys(self):
        'od.iterkeys() -> an iterator over the keys in od'
        return iter(self)

    def itervalues(self):
        'od.itervalues -> an iterator over the values in od'
        for k in self:
            yield self[k]

    def iteritems(self):
        'od.iteritems -> an iterator over the (key, value) items in od'
        for k in self:
            yield (k, self[k])

    def update(*args, **kwds):
        '''od.update(E, **F) -> None.  Update od from dict/iterable E and F.

        If E is a dict instance, does:           for k in E: od[k] = E[k]
        If E has a .keys() method, does:         for k in E.keys(): od[k] = E[k]
        Or if E is an iterable of items, does:   for k, v in E: od[k] = v
        In either case, this is followed by:     for k, v in F.items(): od[k] = v

        '''
        if len(args) > 2:
            raise TypeError('update() takes at most 2 positional '
                            'arguments (%d given)' % (len(args),))
        elif not args:
            raise TypeError('update() takes at least 1 argument (0 given)')
        self = args[0]
        # Make progressively weaker assumptions about "other"
        other = ()
        if len(args) == 2:
            other = args[1]
        if isinstance(other, dict):
            for key in other:
                self[key] = other[key]
        elif hasattr(other, 'keys'):
            for key in other.keys():
                self[key] = other[key]
        else:
            for key, value in other:
                self[key] = value
        for key, value in kwds.items():
            self[key] = value

    __update = update  # let subclasses override update without breaking __init__

    __marker = object()

    def pop(self, key, default=__marker):
        '''od.pop(k[,d]) -> v, remove specified key and return the corresponding value.
        If key is not found, d is returned if given, otherwise KeyError is raised.

        '''
        if key in self:
            result = self[key]
            del self[key]
            return result
        if default is self.__marker:
            raise KeyError(key)
        return default

    def setdefault(self, key, default=None):
        'od.setdefault(k[,d]) -> od.get(k,d), also set od[k]=d if k not in od'
        if key in self:
            return self[key]
        self[key] = default
        return default

    def __repr__(self, _repr_running={}):
        'od.__repr__() <==> repr(od)'
        call_key = id(self), _get_ident()
        if call_key in _repr_running:
            return '...'
        _repr_running[call_key] = 1
        try:
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, self.items())
        finally:
            del _repr_running[call_key]

    def __reduce__(self):
        'Return state information for pickling'
        items = [[k, self[k]] for k in self]
        inst_dict = vars(self).copy()
        for k in vars(OrderedDict()):
            inst_dict.pop(k, None)
        if inst_dict:
            return (self.__class__, (items,), inst_dict)
        return self.__class__, (items,)

    def copy(self):
        'od.copy() -> a shallow copy of od'
        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        '''OD.fromkeys(S[, v]) -> New ordered dictionary with keys from S
        and values equal to v (which defaults to None).

        '''
        d = cls()
        for key in iterable:
            d[key] = value
        return d

    def __eq__(self, other):
        '''od.__eq__(y) <==> od==y.  Comparison to another OD is order-sensitive
        while comparison to a regular mapping is order-insensitive.

        '''
        if isinstance(other, OrderedDict):
            return len(self)==len(other) and self.items() == other.items()
        return dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

    # -- the following methods are only used in Python 2.7 --

    def viewkeys(self):
        "od.viewkeys() -> a set-like object providing a view on od's keys"
        return KeysView(self)

    def viewvalues(self):
        "od.viewvalues() -> an object providing a view on od's values"
        return ValuesView(self)

    def viewitems(self):
        "od.viewitems() -> a set-like object providing a view on od's items"
        return ItemsView(self)
## end of http://code.activestate.com/recipes/576693/ }}}


def roundToNearest(m,x=None):
  def rounding(x):
    x += m/2.0
    mod = x%m
    return x if mod>m*0.98 else x-mod

  return rounding if x is None else rounding(x)

def static_args(f, **static):
  f_args = getargspec(f).args

  @wraps(f)
  def _f(*args):
    args = list(args)
    for fixed,val in static.iteritems():
      args[f_args.index(fixed)] = val
    return f(*args)
  _f.arglist = f_args

  return _f
  
def fcache(f,rounding=None):
  "Caches results of computation f with optional rounding"
  cache = {}

  @wraps(f)
  def _f(*args): 
    if rounding is not None:
      args = str(map(roundToNearest(rounding),args))
      sargs = str(args)
    try:
      if cache.has_key(sargs):
        return cache[sargs]
      else:
        cache[args] = result = f(*args)
        return result
    except TypeError:
      return f(*args)

  return _f

def trace(f):
  "Prints arguments and return values every time f is called"

  argspec = inspect.getargspec(f)

  @wraps(f)
  def _f(*args,**kwargs):
    if args:
      print "Positional args: " + str(args)

    if kwargs:
      print "Keyword args: " + str(kwargs)

    retval = f(*args,**kwargs)
    print "Return value: {}\n".format(retval)
    return retval

  return _f

def grouper(n, iterable, fillvalue=None):
  "Collect data into fixed-length chunks"
  args = [iter(iterable)] * n
  return izip_longest(fillvalue=fillvalue, *args)

def groupat(predicate, iterable, size=2):

  _fill = lambda x: x+(None,)*(size-len(x))
  out = ()
  for x in iterable:
    if predicate(x) and out:
      yield _fill(out)
      out = (x,)
    else:
      out += (x,)

  yield _fill(out)

if __name__=='__main__':
  test_list = 'ABCDEF'
  for x,y in grouper(2,test_list):
    print x,y
