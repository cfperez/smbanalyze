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
