import inspect

from numpy import roots, real
import numpy as np 

try:
  from scipy.optimize import curve_fit
except ImportError:
  from curve_fit import curve_fit

from constants import parameters,kT
from fplot import _subplot
from useful import fix_args, broadcast, static_args
from collections import OrderedDict, Iterable
from operator import or_, itemgetter

class FitError(Exception):
  pass

def fitWLC_masks(x, y, masks, **fitOptions):
  '''
  Return Fit using MMS_rip_region fit function on data in x and y.
  masks: list of boolean arrays which indicate:
    1) the handle (1st element)
    2) the rips (following elements)
  '''
  assert isinstance(x, Iterable)
  assert isinstance(y, Iterable)
  assert map(lambda m: isinstance(m, np.ndarray), masks)
  assert map(lambda m: m.dtype is np.dtype('bool'), masks)

  fitOptions.setdefault('Lc', max(x[masks[-1]]))
  fitOptions.setdefault('fixed', tuple())
  fitOptions['fixed'] += ('K', 'K1', 'Lp1')
  mask_for_mask = combine_masks(masks)
  masks_trimmed = map(lambda m: m[mask_for_mask], masks)
  fit = Fit(x, y, fitfunc=MMS_rip_region_maker(masks_trimmed), 
    mask=mask_for_mask, **fitOptions)
  return fit
  
def fitWLC(x, f, mask=None, **fitOptions):
  "Fit stretching data to WLC model"
  assert len(x) > 0
  assert len(x) == len(f)
  assert mask is None or len(mask) == len(x)
  fitOptions.setdefault('fitfunc', 'MMS')
  fitOptions.setdefault('Lc', max(x)*1.05)
  try:
    if fitOptions['fitfunc'].startswith('MMS'):
      fitOptions.setdefault('fixed', 'K')
  except AttributeError: pass

  fit = Fit(x, f, mask=mask, **fitOptions)
  return fit

############################################################
## Fitting Functions
############################################################
def gauss(x,mu,sigma,A):
  return A*np.exp( -(x-float(mu))**2 / (2*sigma**2) )
gauss.default = {'mu': 1500, 'sigma':500, 'A': 500}

def doublegauss(x,mu,sigma,A,mu2,sigma2,A2):
  return gauss(x,mu,sigma,A)+gauss(x,mu2,sigma2,A2)

@broadcast
def MS(x,Lp,Lc,F0):
  "Marko-Siggia model of worm-like chain"
  x_ = x/float(Lc)
  A = kT(parameters['T'])/Lp
  return A * (0.25/(1-x_)**2 - 0.25 +x_) - F0
MS.default = {'Lp':20.,'Lc':1150.,'F0':0.1}

@broadcast
def MMS(F, Lp, Lc, F0, K):
  "Modified Marko-Siggia model as a function of force"
  f = float(F-F0)* Lp / kT(parameters['T'])
  inverted_roots = roots([1.0, f-0.75, 0.0, -0.25])
  root_index = int(f>=0.75)*2
  root_of_inverted_MS = real(inverted_roots[root_index])
  return Lc * (1 - root_of_inverted_MS + (F-F0)/float(K))
MMS.default = {'Lp':30., 'Lc':1150., 'F0':0.1, 'K': 1200.}
MMS.inverted = True

def MMS_rip(F, Lp, Lc, F0, K, Lp1, Lc1, K1):
  return MMS(F, Lp, Lc, F0, K) + MMS(F, Lp1, Lc1, F0, K1)
MMS_rip.default = dict(MMS.default, Lp1=1.0, Lc1=0., K1=1600.)
MMS_rip.inverted = True

def MMS_rip2(F, Lp, Lc, F0, K, Lp1, Lc1, K1, helix):
  return MMS_rip(F, Lp, Lc-helix, F0, K, Lp1, Lc1, K1)
MMS_rip2.default = dict(MMS_rip.default, helix=2)
MMS_rip.inverted = True

def _is_between(val, range_):
  assert len(range_) == 2
  low, high = range_
  return (val>=low) & (val<=high)

def MMS_rip_maker(handle_limits, upper_limits, split_point=None):
  assert isinstance(handle_limits, Iterable)
  assert isinstance(upper_limits, Iterable)
  assert len(handle_limits) == 2
  assert len(upper_limits) == 2

  if split_point:
    def MMS_rip_global(F, Lp, Lc, F0, K, Lp1, Lc1, K1):
      handle_ext = MMS(F[:split_point], Lp, Lc, F0, K)
      upper_ext = MMS(F[split_point:], Lp, Lc, F0, K) + \
                  MMS(F[split_point:], Lp1, Lc1, F0, K1)
      return np.append(handle_ext, upper_ext)
  else:      
    def MMS_rip_global(F, Lp, Lc, F0, K, Lp1, Lc1, K1):
      handle_ext = MMS(F, Lp, Lc, F0, K)
      upper_ext = handle_ext + MMS(F, Lp1, Lc1, F0, K1)
      return np.where(_is_between(F, upper_limits), upper_ext, handle_ext)
  MMS_rip_global.default = dict(MMS_rip.default)
  MMS_rip_global.inverted = True

  return MMS_rip_global

def fitWLCrip(x_handle, f_handle, x_upper, f_upper, mask=None, **fitOptions):
  fitOptions.setdefault('Lc', max(x_upper))
  fitOptions.setdefault('fixed', tuple())
  fitOptions['fixed'] += ('K', 'K1', 'Lp1')
  ext = np.append(x_handle, x_upper)
  force = np.append(f_handle, f_upper)
  def limits(data):
    return (min(data), max(data))
  fit = Fit(ext, force, fitfunc=MMS_rip_maker(limits(f_handle), limits(f_upper), split_point=len(x_handle)),
            mask=mask, 
            **fitOptions)
  return fit

def combine_masks(masks):
  if len(masks) == 1:
    return masks[0]
  return np.array(reduce(or_, masks))



def MMS_rip_region_maker(masks):
  '''Creates a fitting function for an arbitrary number of fit regions in masks
  that allows simultaneous fitting of each as MMS rips.
  '''
  assert isinstance(masks, (tuple,list))
  assert len(masks) > 1
  assert map(lambda e: isinstance(e, np.ndarray), masks)

  handle = masks.pop(0)

  def MMS_rip_region(F, Lp, Lc, F0, K, Lp1, Lc1, K1, **rips):
    rips['Lc1'] = Lc1
    handle_force = F[handle]
    rip_forces = [F[mask] for mask in masks]
    rip_sizes = map(lambda r: rips[r], sorted(rips))
    rip_items = zip(rip_forces, rip_sizes)

    handle_ext = MMS(handle_force, Lp, Lc, F0, K)
    rip_ext = [MMS_rip(force, Lp, Lc, F0, K, Lp1, Lc_rip, K1) 
                for force,Lc_rip in rip_items]

    if len(rip_ext) > 1:
      rip_ext = np.append(*rip_ext)
    return np.append(handle_ext, rip_ext)
  
  addl_rips = ['Lc{}'.format(n) for n in range(2,1+len(masks))]
  MMS_rip_region.default = MMS_rip.default.copy()
  MMS_rip_region.default.update([(rip,10) for rip in addl_rips])
  MMS_rip_region.arglist = ['F', 'Lp', 'Lc', 'F0', 'K', 'Lp1', 'Lc1', 'K1'] + addl_rips
  MMS_rip_region.inverted = True

  return MMS_rip_region


class Fit(object):
  @classmethod
  def extends(cls, fit, x, y, fitfunc, fixed=(), mask=None, verbose=False, **user_parameters):
    fixed = tuple(fixed) + tuple(fit.parameters)
    params = fit.parameters.copy()
    params.update(user_parameters)
    newfit = cls(x, y, fitfunc, fixed, mask, verbose=verbose, **params)
    return newfit
    
  def __init__(self, x, y, fitfunc, fixed=(), mask=None, weights=None, verbose=False, **user_parameters):
    "Initialize to a specific fitting function, optionally fitting to data specified"

    if isinstance(fitfunc,str):
      fitfunc = eval(fitfunc)
    self.fitfunc = fitfunc
    self.inverted = getattr(fitfunc, 'inverted', False)

    if mask is not None:
      self.mask = np.logical_not(mask)
      to_masked = lambda ar: ar[mask]
    else:
      self.mask = None
      to_masked = lambda ar: ar

    # Use inspection to get parameter names from fit function
    # assuming first argument is independent variable
    arg_names = getattr(fitfunc, 'arglist', inspect.getargspec(fitfunc).args)[1:]
    fit_parameters = OrderedDict.fromkeys(arg_names)
    valid_args = frozenset(arg_names)

    if not set(user_parameters.keys()) <= valid_args:
      raise FitError('Keyword arguments {0} can only set valid fit parameters {1}'.format(user_parameters.keys(), list(valid_args)))

    # Don't bother getting function defaults of user specifies everything
    if set(user_parameters.keys()) == valid_args:
      fit_parameters.update(user_parameters)
    else:
      try:
        fit_defaults = getattr(fitfunc, 'default', {}).items() or \
          zip( reversed(arg_names), reversed(fitfunc.func_defaults))
        fit_parameters.update(fit_defaults, **user_parameters)
      except TypeError:
        raise FitError("Missing function defaults. Must specify in **user_parameters")

    if isinstance(fixed, str):
        fixed = (fixed,)
    if not set(fixed) <= valid_args:
      raise FitError('Fixed argument(s) [{}] must specify one of: {}'.format(
            fixed, valid_args))

    free_parameters = valid_args - set(fixed)

    starting_p = [ param for key,param in fit_parameters.iteritems()
        if key not in fixed ]

    if fixed:
      fixed_parameters = OrderedDict( filter(lambda item: item[0] in fixed,
                                fit_parameters.items()) )
      fitfunc = fix_args(fitfunc, **fixed_parameters)

    x, y = to_masked(x), to_masked(y)
    if self.inverted:
      x,y = y,x
    self.x = x

    param_best_fit, self.covariance = curve_fit(fitfunc, x, y, starting_p, sigma=weights)
    self.fitOutput = fitfunc(x, *param_best_fit)
    self.residual = self.fitOutput-y

    self.fixed = fixed
    self.parameters = fit_parameters.copy()
    self.free_parameters = OrderedDict({})

    try:
      param_best_fit = param_best_fit.tolist()
      for param in fit_parameters:
        if param in free_parameters:
          v = param_best_fit.pop(0)
          self.free_parameters[param] = v
          self.parameters[param] = v
    except IndexError:
      raise FitError("Free/fix parameter mismatch!")
    error = self.covariance.diagonal()
    self.error = OrderedDict(zip(self.free_parameters.keys(), error))

  def __call__(self, x=None):
    return self.fitfunc(x, *self.parameters.values())
	
  def __getitem__(self,key):
    return self.parameters[key]

  def __len__(self):
    return len(self.parameters)

  def _to_plot(self, **kwargs):
    x = sorted(self.x)
    y = self(x)
    if self.inverted:
      x,y = y,x
    kwargs.setdefault('linewidth', 2)
    return (x,y), kwargs

  def plot(self, **kwargs):
    args = (self.x, self.fitOutput)
    if self.inverted:
      args = reversed(args)
    kwargs.setdefault('linewidth', 2)
    return _subplot(*args,**kwargs)

  def __repr__(self):
    return "<Fit function '{0}' using parameters {1}".format(self.fitfunc.func_name, 
        ', '.join(['{0}={1:.2f}'.format(*p) for p in self.parameters.items()]))

  def __unicode__(self):
    params = self.free_parameters
    return self.fitfunc.func_name + ' fit: ' + \
      ' '.join(u'{0}={1:.2f}\u00B1{2:.2f}'.format(p,params[p],self.error[p]) 
                                            for p in params)

  def __str__(self):
    return unicode(self).encode('utf-8')

  def toFile(self,filename):
    raise NotImplementedError

class FitFret(Fit):
  fitfunc = gauss
  sigma = 0.2

  def __init__(self, x, y, numGauss):
    self.fitfunc = funcBuilder( gauss, numGauss)
    self.fitfunc.params = ('mu','sigma','A')*numGauss

    guess = np.linspace(0.1,0.9,numGauss)
    A = max(y)
    self.x = np.linspace(0,1)

    p0 = []
    for mu in guess:
      p0.extend([mu,FitFret.sigma,A])

    super(FitFret,self).__init__(self.fitfunc,x,y,p0)

  def plot(self,**kwargs):
    return super(FitFret,self).plot(self.x,**kwargs)


def funcBuilder(func,num):
  numArgs = len(inspect.getargspec(func).args)-1

  def g(*args):
    x=args[0]
    out = 0
    for i in range(1,num*numArgs,numArgs):
      out += func(x,*args[i:i+numArgs])
    return out

  return g

