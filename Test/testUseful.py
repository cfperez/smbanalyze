import inspect
from smbanalyze.useful import fix_args
from smbanalyze.curvefit import MS, MMS, MMS_rip

def testReturnsFunction():
  f = lambda x: x
  assert inspect.isfunction( fix_args(f) )

def testReturnsIdentityCase():
  f = lambda x: x
  result = fix_args(f)
  assert result(100) == f(100)

def testOneConstantArg():
  f = lambda x,y: x+y
  result = fix_args(f,y=10)
  assert result(100) == f(100,10)

def testTwoConstantArgs():
  f = lambda x,y,z: (x,y,z)
  result = fix_args(f,y=50,z=50)
  assert result(100) == f(100,50,50)

def testTwoConstantArgsSkipFirst():
  f = lambda x,y,z: (x,y,z)
  result = fix_args(f,x=50,z=50)
  assert result(100) == f(50,100,50)

def testTwoConstantArgsWithKeywords():
  def f(x, p1, p2=10):
    return x, p1, p2
  result = fix_args(f, p1=20, p2=30)
  assert result(66) == f(66, 20, 30)

def testMSfixedparams():
  for x in range(100,1000,100):
    yield MSfitfunc, x

def MSfitfunc(X):
  to_fix = {'Lp': 30, 'Lc': 1100, 'F0': 0}
  fixed = fix_args(MS, **to_fix)
  assert fixed(X) == MS(X, **to_fix)

def testMMSfixedParams():
  for x in range(100,1000,100):
    yield MMSfitfunc, x
  
def MMSfitfunc(X):
  to_fix = {'Lp': 30, 'Lc': 1100, 'F0': 0, 'K': 1200}
  fixed = fix_args(MMS, **to_fix)
  assert fixed(X) == MMS(X, **to_fix)

if __name__=="__main__":
  unittest.main()
