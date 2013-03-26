import inspect
from operator import isSequenceType
from itertools import izip

from numpy import arccos, cos, exp, sqrt, log, fabs, pi, roots, real
import numpy as np 

try:
  from scipy.optimize import curve_fit
except ImportError:
  from curve_fit import curve_fit

from constants import parameters,kT
from fplot import _subplot
from useful import OrderedDict, fix_args, broadcast

class FitError(Exception):
  pass

def fitWLC(x, f, mask=None, **fitOptions):
  "Fit stretching data to WLC model"
  fitOptions.setdefault('fitfunc', 'MMS')
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

def doublegauss(x,mu,sigma,A,mu2,sigma2,A2):
  return gauss(x,mu,sigma,A)+gauss(x,mu2,sigma2,A2)

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
  inverted_roots = roots([1, f-0.75, 0, -0.25])
  root_index = int(f>=0.75)*2
  root_of_inverted_MS = real(inverted_roots[root_index])
  return Lc * (1 - root_of_inverted_MS + (F-F0)/float(K))
MMS.default = {'Lp':30., 'Lc':1150., 'F0':0.1, 'K': 1200.}
MMS.inverted = True

def MMS_rip(F, Lp=30., Lc=1150., F0=0.1, K=1200., Lp1=2.5, Lc1=0., K1=1100.):
  return MMS(F, Lp, Lc, F0, K) + MMS(F, Lp1, Lc1, F0, K1)
MMS_rip.default = dict(MMS.default, Lp1=2.5, Lc1=0., K1=1100.)
MMS_rip.inverted = True

############################################################
## Fit class
############################################################
class Fit(object):
  def __init__(self, x, y, fitfunc, fixed=(), mask=None, verbose=False, **user_parameters):
    "Initialize to a specific fitting function, optionally fitting to data specified"

    if isinstance(fitfunc,str):
      fitfunc = eval(fitfunc)
    self.fitfunc = fitfunc
    self.inverted = getattr(fitfunc, 'inverted', False)

    self.mask = np.logical_not(mask)
    #to_masked = lambda ar: np.ma.array(ar, self.mask)
    to_masked = lambda ar: ar[mask]

    # Use inspection to get parameter names from fit function
    # assuming first argument is independent variable
    arg_names = getattr(fitfunc, 'arglist', inspect.getargspec(fitfunc).args)[1:]
    fit_parameters = OrderedDict.fromkeys(arg_names)
    valid_args = frozenset(arg_names)

    if not set(user_parameters.keys()) <= valid_args:
      raise FitError('Keyword arguments can only set valid fit parameters')

    fit_defaults = getattr(fitfunc, 'default', {}).items() or \
      zip( reversed(arg_names), reversed(fitfunc.func_defaults))
    fit_parameters.update(fit_defaults, **user_parameters)

    if isinstance(fixed, str):
      fixed = (fixed,)
    elif not isSequenceType(fixed):
      raise ValueError("Argument 'fixed' must be a string or a tuple: %s" % fixed)
    elif not set(fixed) <= valid_args:
      raise FitError('Fixed argument must specify one of: %s' % valid_args)

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

    fit_params, self.error = curve_fit(fitfunc, x, y, starting_p)
    self.fitOutput = fitfunc(x, *fit_params)
    self.residual = self.fitOutput-y

    self.fixed = fixed
    self.parameters = fit_parameters.copy()
    try:
      fit_params = fit_params.tolist()
      for param in fit_parameters:
        if param in free_parameters:
          self.parameters[param] = fit_params.pop(0)
    except IndexError:
      raise FitError("Free/fix parameter mismatch!")

  def __call__(self, x=None):
    return self.fitfunc(x, *self.parameters.values())
	
  def __getitem__(self,key):
    return self.parameters[key]

  def __len__(self):
    return len(self.parameters)

  def _to_plot(self, **kwargs):
    if self.inverted:
      x,y = self.fitOutput,self.x
    else:
      x,y = self.x,self.fitOutput
    kwargs.setdefault('linewidth', 2)
    return (x,y), kwargs

  def plot(self, **kwargs):
    args, kwargs = self._to_plot(**kwargs)
    return _subplot(*args,**kwargs)

  def __repr__(self):
    return "<Fit function '{0}' using parameters {1}".format(self.fitfunc.func_name, 
        ', '.join(['{0}={1:.2f}'.format(*p) for p in self.parameters.items()]))

  def __str__(self):
    return self.fitfunc.func_name + ' fit: ' + \
      ' '.join('{0}={1:.2f}'.format(p,v) for p,v in self.parameters.items())

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

