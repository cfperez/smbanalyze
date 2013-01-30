import inspect
from functools import wraps
from itertools import izip
from numpy import arccos, cos, exp, sqrt, log, fabs, pi

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
  "MS(x, Lp, Lc, K) = Marko-Siggia model of worm-like chain"
  x_ = x/float(Lc)
  A = kT(parameters['T'])/Lp
  return A * (0.25/(1-x_)**2 - 0.25 +x_) + F0
MS.default = {'Lp':20,'Lc':1100,'F0':0.1}

def MMS2(F,Lp,Lc,F0,K):
  "Calculates x(F) for Modified Marko-Siggia"

  #if not np.iterable(F):
  #  F = [F]
  #F = np.asarray(F)

  # Normalize parameters
  Lp_ = Lp/kT(parameters['T'])
  F_ = (F-F0)/K

  # Create variables for computation
  Q = 1 + F_
  P = F_ + Lp_*F + 0.25

  #U = -(Lp_*F - 0.75)**2 / 9.
  #V = -(Lp_*F - 0.75)**3 / 27. + 0.125
  #L_ = -(2*Q + P)/3.
  L = -(2*Q + P)/3.

  u = -(P-Q)**2 / 9.0   #=(-1/9)*(Lp_*F - 3/4)**2
  v = -(P-Q)**3 / 27.0 + 0.125   #=(-1/27)*(Lp_*F - 3/4)**3 + 1/8

  #def calc(u,v,L):
  if v**2+u**3 < 0.0:
    theta = arccos(sqrt(-v**2/u**3)) / 3.0
    if v<0:
      return (2*sqrt(-u)*cos(theta+2*pi/3.) - L)*Lc
    else:
      return (-2*sqrt(-u)*cos(theta) - L)*Lc
  else:
    A = exp( log(-v+sqrt(0.015625 - ((Lp_*F-0.75)**3/108.) ) ) )
    B = exp( log(-v-sqrt(0.015625 - ((Lp_*F-0.75)**3/108.) ) ) )
    return -Lc*(fabs(A+B) + L)

  #return np.array([calc(u,v,L) for u,v,L in izip(U,V,L_)])
  #return calc(U,V,L_)
MMS2.default = {'Lp':20,'Lc':1100,'F0':0.1, 'K': 1100}
MMS2.inverted = True

def MMS(x,Lp,Lc,K,startf=20.):
  x_ = x/float(Lc)
  guess = [startf/K]

  A = kT(parameters['T'])/Lp
  wlc = lambda F_: A*(0.25/(1-x_+F_)**2 - 0.25 + x_ - F_)
  wlc_min = lambda F_: (np.asarray(F_)*K - wlc(F_))**2
  wlc_deriv = lambda F_: 2*(F_*K-wlc(F_)) * (K - A*(0.5/(1-x_-F_)**3-1))
  return fmin_bfgs(wlc_min,guess,fprime=wlc_deriv,maxiter=200,disp=False)[0]*K
MMS.default = MS.default

def minimize(func,guess,*args,**kw):
  return np.array([func(x) for x in guess])
  #return np.asarray([fmin_bfgs(func,[x],disp=False,*args,**kw)[0] for x in guess])

def funcBuilder(func,num):
  numArgs = len(inspect.getargspec(func).args)-1

  def g(*args):
    x=args[0]
    out = 0
    for i in range(1,num*numArgs,numArgs):
      out += func(x,*args[i:i+numArgs])
    return out

  return g

def apply(fitfunc,*args):
  return lambda x: fitfunc(x,*args)

def broadcast(f,var=None):

  @wraps(f)
  def _f(X,*args):
    if np.iterable(X):
      return np.array([f(x,*args) for x in X])
    else:
      return f(X,*args)
  _f.arglist = inspect.getargspec(f).args

  return _f

############################################################
## Fit class
############################################################
class Fit(object):
  def __init__(self, fitfunc, x, y, **kwargs):
    "Initialize to a specific fitting function, optionally fitting to data specified"

    verbose = kwargs.get('verbose',True)

    self.x = x

    if isinstance(fitfunc,str):
      fitfunc = eval(fitfunc)

    # Use inspection to get parameter names from fit function
    # assuming first argument is independent variable
    try:
      arg_names = getattr(fitfunc,'arglist',inspect.getargspec(fitfunc).args)[1:]
    except TypeError:
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
      to_fix = dict( (param,kwargs[param]) for param in fixed )
      self.__dict__.update(to_fix)
      self.fitfunc = useful.fix_args(fitfunc, **to_fix)
    else:
      self.fitfunc = fitfunc

    if getattr(fitfunc,'inverted',False):
      x,y=y,x

    fit_params,self.error = curve_fit(self.fitfunc,x,y,starting_param)
    self.fitY = self.fitfunc(x,*fit_params)
    self.residual = self.fitY-y

    self.free_parameters = dict( zip(free_params,fit_params) )
    self.__dict__.update(self.free_parameters)

    self.parameters = dict(self.free_parameters).update(to_fix)

    if verbose:
      print "Fitting with parameters {0}".format(','.join(['{0}={1:.2f}'.format(*p) for p
        in zip(free_params,fit_params)]))

  def __call__(self,x=None):
	return self.fitfunc(x,*self.free_parameters.values)
	
  def __getitem__(self,key):
    return self.param_names.index(key)

  def __len__(self):
	return len(self.param_names)

  def plot(self,x=None,**kwargs):
	if x is None:
	  x = self.x
	return plt.plot(x,self(x),**kwargs)

  def __repr__(self):
	if hasattr(self.fitfunc,'params'):
	  return self.fitfunc.func_name + ': ' + \
		'; '.join(name+'='+str(value) for name,value in zip(self.fitfunc.params,self.params))
	else:
	  return '<Fit {0} using params {1}'.format(self.fitfunc.func_name,list(self))

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
