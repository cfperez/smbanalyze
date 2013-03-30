import numpy as np

from smbanalyze import curvefit

X = None
Y = None
params = (27, 1046, 0, 1200)

def setUp():
  pass
  global X,F,fit_params
  X = range(700,1100,10)
  F = range(0,21,1)
  fit_params = (30, 1150, 0, 1e6)
  #data = FileIO.loadPull('/Users/blocklab/Analysis/2012.09.06/SJF4_1B_0.5nM_s1m6.str')
  #X = data.ext[:152]
  #Y = data.f[:152]
  
def testMSandMMSequivalent():
  msForce = curvefit.MS(X, *fit_params[:-1])
  mmsExt = curvefit.MMS(msForce, *fit_params)
  assert max(abs(mmsExt-X)) < 0.05

def testMMSreturntype():
  pass
  #for x in X[0::20]:
  #  yield checkMMSreturntype, x

def checkMMSreturntype(x):
  pass
  #eval = curvefit.MMS(x,*params)
  #assert isinstance(eval,type(x))

def testMMSreturnvalue():
  pass
  #for x,y in zip(X[::20],Y[::20]):
  #  yield checkMMSeval,x,y

def testMMSvectorinput():
  pass
  #for x,y in zip(X[::10],Y[::10]):
    #yield checkMMSeval,[x,x],[y,y]

def checkMMSeval(x,y):
  pass
  #results = curvefit.MMS(x,*params)
  #assert np.all(abs(results-y < 10))

def testMSfitHasParams():
  pass
  #fit = curvefit.Fit(curvefit.MS,X,Y)
  #for p in curvefit.MS.default:
    #assert hasattr(fit,p)

def testMMSfitHasParams():
  pass
  #fit = curvefit.Fit(curvefit.MMS,X,Y)
  #for p in curvefit.MMS.default:
    #assert hasattr(fit,p)
  
