import numpy as np

import curvefit
import FileIO

X = None
Y = None
params = [27, 1046, 1200]

def setUp():
  global X,Y
  data = FileIO.loadPull('/Users/blocklab/Analysis/2012.09.06/SJF4_1B_0.5nM_s1m6.str')
  X = data.ext[:152]
  Y = data.f[:152]
  
def testMMSreturntype():
  for x in X[0::20]:
    yield checkMMSreturntype, x

def checkMMSreturntype(x):
  eval = curvefit.MMS(x,*params)
  assert isinstance(eval,type(x))

def testMMSreturnvalue():
  for x,y in zip(X[::20],Y[::20]):
    yield checkMMSeval,x,y

def testMMSvectorinput():
  for x,y in zip(X[::10],Y[::10]):
    yield checkMMSeval,[x,x],[y,y]

def checkMMSeval(x,y):
  results = curvefit.MMS(x,*params)
  assert np.all(abs(results-y < 10))

def testMSfitHasParams():
  fit = curvefit.Fit(curvefit.MS,X,Y)
  for p in curvefit.MS.default:
    assert hasattr(fit,p)

def testMMSfitHasParams():
  fit = curvefit.Fit(curvefit.MMS,X,Y)
  for p in curvefit.MMS.default:
    assert hasattr(fit,p)
  
