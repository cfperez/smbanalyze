from scipy.optimize import curve_fit
import numpy as np 
import matplotlib.pyplot as plt
import inspect

def gauss(x,mu,sigma,A):
  return A*np.exp( -(x-float(mu))**2 / (2*sigma**2) )
gauss.params = ('mu','sigma','A')

def doublegauss(x,mu,sigma,A,mu2,sigma2,A2):
  return gauss(x,mu,sigma,A)+gauss(x,mu2,sigma2,A2)

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

def fit( fitfunc, x, y, guess):
  p0 = [guess,0.2,10]
  return Fit(fitfunc.func_name, *curve_fit( fitfunc, x, y, p0=p0 ) )

class Fit(object):
  def __init__(self, fitfunc, x, y, p0):
	self.fitfunc = fitfunc
	self.params,self.error = curve_fit(fitfunc,x,y,p0)
	for i,p in enumerate(fitfunc.params):
	  setattr(self,p,self.params[i])

  def __call__(self,x):
	return self.fitfunc(x,*self.params)
	
  def __getitem__(self,key):
	return self.params[key]

  def __len__(self):
	return len(self.params)

  def plot(self,x=None,**kwargs):
	if not x:
	  x = np.linspace(0,1)

	return plt.plot(x,self(x),**kwargs)

  def __repr__(self):
	return self.fitfunc.func_name + ': ' + '; '.join(name+'='+str(value) for name,value in zip(self.fitfunc.params,self.params))

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
