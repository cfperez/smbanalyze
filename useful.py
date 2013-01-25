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
  return int(s) if s else None

def isInt(s):
  try:
    return s == str(toInt(s))
  except ValueError:
    return False

def negate(it):
  return [-i for i in it]

class dotdict(dict):
  def __init__(self,*args,**kwargs):
	super(dotdict,self).__init__(*args,**kwargs)
	self._locked = False

  def __setitem__(self, key, val):
	if not str(key).replace('_','').isalnum():
	  raise KeyError, "Key must be alphanumeric"
	super(dotdict,self).__setitem__(key,val)

  def __getitem__(self,key):
	if not self.has_key(key) and not str(key).startswith('_') and not self._locked:
	  self[key] = dotdict()
	return super(dotdict,self).__getitem__(key)

  def __getattr__(self, name):
	if not self.__dict__.has_key(name):
	  try:
		return self.__getitem__(name)
	  except KeyError:
		raise AttributeError('Dotdict has no attribute %s' % name)
	else: return super(dotdict,self).__getattr__(name)

  def __setattr__(self, name, value):
	if name.startswith('_'):
	  super(dotdict,self).__setattr__(name,value)
	elif not self._locked:
	  self[name] = value
	return self

  def __iter__(self):
	for val in self.values():
	  yield val

  def __iadd__(self,other):
    for key in other:
      if self.has_key(key):
        self[key] = np.append(self[key],other[key])
      else:
        self[key] = other[key]
      return self

  def __add__(self,other):
	temp = dotdict(self)
	temp.update(other)
	return temp

  def _lock(self):
    for value in self:
      if hasattr(value,'_lock'):
        value._lock()
    self._locked = True

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
  
def fix_args(f,**fixed):
  "Return function f "
  f_args = getargspec(f).args
  fixed_i = map(f_args.index,fixed.keys())

  @wraps(f)
  def _f(*args):
    args = list(args)
    for i,fixed_arg in zip(fixed_i,fixed.values()):
      args.insert(i,fixed_arg)
    return f(*args)
  _f.arglist=f_args

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
