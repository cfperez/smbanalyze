import inspect

from numpy import roots, real
import numpy as np 

try:
  from scipy.optimize import curve_fit
except ImportError:
  from curve_fit import curve_fit

from constants import parameters,kT
from fplot import _subplot
from useful import fix_args, broadcast
from collections import OrderedDict

class FitError(Exception):
  pass

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

############################################################
## Fit class
############################################################
class Fit(object):
  @classmethod
  def extends(cls, fit, x, y, fitfunc, fixed=(), mask=None, verbose=False, **user_parameters):
    fixed = tuple(fixed) + tuple(fit.parameters)
    params = fit.parameters.copy()
    params.update(user_parameters)
    newfit = cls(x, y, fitfunc, fixed, mask, verbose, **params)
    #newfit.error.update(fit.error)
    return newfit
    
  def __init__(self, x, y, fitfunc, fixed=(), mask=None, verbose=False, **user_parameters):
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

    param_best_fit, self.covariance = curve_fit(fitfunc, x, y, starting_p)
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
    args, kwargs = self._to_plot(**kwargs)
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

