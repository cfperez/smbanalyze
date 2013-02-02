import inspect
from functools import wraps, partial, update_wrapper
from itertools import izip
from numpy import arccos, cos, exp, sqrt, log, fabs, pi, roots, real

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,fmin_bfgs
import numpy as np 

from Constants import parameters,kT
import useful

class FitError(Exception):
  pass

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
  return A * (0.25/(1-x_)**2 - 0.25 +x_) + F0
MS.default = {'Lp':20.,'Lc':1150.,'F0':0.1}

@useful.broadcast
def MSinv(F, Lp=20.0, Lc=1150.0, F0=0.1, K=1200.0):
  f = float(F-F0)* Lp / kT(parameters['T'])
  inverted_roots = roots([1, f-0.75, 0, -0.25])
  root_index = int(f>=0.75)*2
  root_of_inverted_MS = real(inverted_roots[root_index])
  return Lc * (1 - root_of_inverted_MS + (F-F0)/float(K))
MSinv.default = {'Lp':20., 'Lc':1150., 'F0':0.1, 'K': 1200.}
MSinv.inverted = True
MSinv.arglist = ('F','Lp','Lc','F0','K')

############################################################
## Fit class
############################################################
class Fit(object):
  def __init__(self, fitfunc, x, y, **kwargs):
    "Initialize to a specific fitting function, optionally fitting to data specified"

    verbose = kwargs.get('verbose',True)

    if isinstance(fitfunc,str):
      fitfunc = eval(fitfunc)

    # Use inspection to get parameter names from fit function
    # assuming first argument is independent variable
    try:
      arg_names = getattr(fitfunc,'arglist',inspect.getargspec(fitfunc).args)[1:]
    except AttributeError:
      arg_names = kwargs.pop('parameters')
    fixed = kwargs.pop('fixed',())
    if not isinstance(fixed,tuple):
      fixed = (fixed,)

    # Parameter names (from inspection) that we will fit to (not fixed)
    free_params = [ p for p in arg_names if p not in fixed ]

    # Use parameter defaults defined by the fitfunc itself
    default = getattr(fitfunc,'default', {})

    starting_param = [ kwargs.setdefault(p,default.get(p,1)) for p in free_params ]

    if fixed:
      to_fix = useful.OrderedDict( (param,kwargs.get(param,default[param])) for param in fixed )
      self.__dict__.update(to_fix)
      self.fitfunc = partial(fitfunc, **to_fix)
    else:
      self.fitfunc = fitfunc
      to_fix = {}

    self.inverted = getattr(fitfunc, 'inverted', False)
    if self.inverted:
      x,y=y,x
    self.x = x

    fit_params,self.error = curve_fit(self.fitfunc,x,y,starting_param)
    self.fitOutput = self.fitfunc(x,*fit_params)
    self.residual = self.fitOutput-y

    self.free_parameters = useful.OrderedDict( zip(free_params,fit_params) )
    self.__dict__.update(self.free_parameters)

    self.parameters = useful.OrderedDict(self.free_parameters)
    if to_fix:
      self.parameters.update(to_fix)

    if verbose:
      print "Fitting with parameters {0}".format(','.join(['{0}={1:.2f}'.format(*p) for p
        in zip(self.free_parameters,fit_params)]))

  def __call__(self,x=None):
    return self.fitfunc(x,*self.free_parameters.values())
	
  def __getitem__(self,key):
    return self.parameters[key]

  def __len__(self):
    return len(self.parameters)

  def plot(self,**kwargs):
    if self.inverted:
      x,y = self.fitOutput,self.x
    else:
      x,y = self.x,self.fitOutput
    return plt.plot(x,y,**kwargs)

  def __repr__(self):
    return '<Fit {0} using parameters {1}'.format(
      self.fitfunc.func_name, self.parameters)

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

