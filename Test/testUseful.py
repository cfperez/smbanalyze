import unittest
import inspect
from useful import fix_args

class TestFixArgs(unittest.TestCase):
  
  def testReturnsFunction(self):
    f = lambda x: x
    assert inspect.isfunction( fix_args(f) )

  def testReturnsIdentityCase(self):
    f = lambda x: x
    result = fix_args(f)
    assert result(100) == f(100)

  def testOneConstantArg(self):
    f = lambda x,y: x+y
    result = fix_args(f,y=10)
    assert result(100) == f(100,10)

  def testTwoConstantArgs(self):
    f = lambda x,y,z: (x,y,z)
    result = fix_args(f,x=50,z=50)
    assert result(100) == f(50,100,50)

if __name__=="__main__":
  unittest.main()
