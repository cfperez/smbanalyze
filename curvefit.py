import inspect
import numpy as np 

try:
  from scipy.optimize import curve_fit
except ImportError:
  from curve_fit import curve_fit

import progressbar as pbar
from fplot import subplot
from useful import fix_args
from collections import OrderedDict, Iterable
from operator import or_

class FitError(Exception):
  pass

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
  
def convert_masks_to_contiguous_regions(masks, mask_for_mask):
  return map(lambda m: m[mask_for_mask], masks)

def combine_masks_with_or(masks):
  if len(masks) == 1:
    return masks[0]
  return np.array(reduce(or_, masks))

############################################################
## Fitting Functions
############################################################
def gauss(x,mu,sigma,A):
  return A*np.exp( -(x-float(mu))**2 / (2*sigma**2) )
gauss.default = {'mu': 1500, 'sigma':500, 'A': 500}

def doublegauss(x,mu,sigma,A,mu2,sigma2,A2):
  return gauss(x,mu,sigma,A)+gauss(x,mu2,sigma2,A2)

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


class Fit(object):
  @classmethod
  def extends(cls, fit, x, y, fitfunc, fixed=(), mask=None, **user_parameters):
    fixed = tuple(fixed) + tuple(fit.parameters)
    params = fit.parameters.copy()
    params.update(user_parameters)
    newfit = cls(x, y, fitfunc, fixed, mask, **params)
    return newfit
    
  def __init__(self, x, y, fitfunc, fixed=(), mask=None, weights=None, maxfev=None,
    **user_parameters):
    "Initialize to a specific fitting function, optionally fitting to data specified"

    if isinstance(fitfunc,str):
      fitfunc = eval(fitfunc)
    self.fitfunc = fitfunc
    self.inverted = getattr(fitfunc, 'inverted', False)

    if mask is not None:
      self.mask = mask
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

    # Don't bother getting function defaults if user specifies everything
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

    maxfeval = maxfev or int(np.sqrt(len(x)/3)+50)
    with pbar.done_on_complete(maxfeval, status='Fitting:') as auto_pbar:
      param_best_fit, self.covariance = curve_fit(
        pbar.progress_on_call(fitfunc, auto_pbar), x, y, starting_p, sigma=weights, maxfev=maxfeval)
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

  def __call__(self, x):
    return self.fitfunc(x, *self.parameters.values())
	
  def __iter__(self):
    return iter(self.parameters)

  def __getitem__(self, key):
    return self.parameters[key]

  @property
  def reduced_chisquared(self):
    DOF = len(self.fitOutput)-len(self.free_parameters)-1
    return sum(self.residual**2)/np.var(self.residual)/DOF

  def _to_plot(self, **kwargs):
    x = sorted(self.x)
    y = self(x)
    if self.inverted:
      x,y = y,x
    return (x,y), kwargs

  def plot(self, style='.', **kwargs):
    x,y = (self.x,self.fitOutput) if not self.inverted else (self.fitOutput,self.x)
    #kwargs.setdefault('marker', '.')
    kwargs.setdefault('markersize', 3)
    #kwargs.setdefault('linestyle', '')
    return subplot(x, y, style, **kwargs)

  def __getstate__(self):
    state = self.__dict__.copy()
    state['fitfunc'] = (self.fitfunc.__module__, self.fitfunc.func_name)
    return state

  def __setstate__(self, state):
    self.__dict__ = state
    func_module, func_name = state['fitfunc']
    import importlib
    func_mod = importlib.import_module(func_module)
    self.fitfunc = getattr(func_mod, func_name)

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

class FitRegions(Fit):
  def __init__(self, x, y, fitfunc, fixed=(), regions=None, weights=None, **user_parameters):
    self.fitfunc_generator = fitfunc
    fitfunc, combined_mask = self._make_fitfunc_from_regions(fitfunc, regions)
    super(FitRegions, self).__init__(
        x, y, 
        fitfunc,
        mask=combined_mask,
        fixed=fixed,
        weights=weights,
        **user_parameters
    )
    self.regions = regions

  def _make_fitfunc_from_regions(self, func_gen, regions):
    combined_mask = combine_masks_with_or(regions)
    masks_trimmed = convert_masks_to_contiguous_regions(regions, combined_mask)
    return func_gen(masks_trimmed), combined_mask

  def __setstate__(self, state):
    self.__dict__ = state
    fitfunc_maker = state['fitfunc_generator']
    self.fitfunc, mask = self._make_fitfunc_from_regions(fitfunc_maker, self.regions)


