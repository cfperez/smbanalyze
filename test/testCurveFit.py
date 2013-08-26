import numpy as np
from smbanalyze import curvefit
import unittest

X = None
Y = None
params = (27, 1046, 0, 1200)

class fitWLCripTest(unittest.TestCase):

    DEFAULT_HANDLE_LIMIT = (1,10)
    DEFAULT_UPPER_LIMIT = (13,20)

    def setUp(self):
        self.handle_params = (30, 1150, 0, 1200)
        self.rip_test_ext, self.rip_test_force = self._generate_test_data()

    def tearDown(self):
        pass

    def rip_params(self, Lc=10):
        return self.handle_params + (1, Lc, 1600)

    def _generate_test_data(self, rip_sizes=[10],
                            handle_limit=DEFAULT_HANDLE_LIMIT, 
                            upper_limit=DEFAULT_UPPER_LIMIT):
        f_handle, f_upper = np.linspace(*handle_limit), np.linspace(*upper_limit)
        ext, forces = self._generate_test_regions([handle_limit, upper_limit], rip_sizes)
        return np.append(*ext), np.append(*forces)

    def _generate_test_regions(self, limits, Lc=[]):
        self.assertGreaterEqual(len(limits), 1)
        self.assertGreaterEqual(len(limits), len(Lc)+1,
            msg='Must have limits for handle and every rip specified by arg "Lc"')
        forces = map(lambda x: np.linspace(*x), limits)
        ext_out = [curvefit.MMS(forces[0], *self.handle_params)]
        if Lc:
            ext_out.append(*[curvefit.MMS_rip(force, *self.rip_params(l)) 
                for force,l in zip(forces[1:], Lc)])
        return ext_out, forces

    def _generate_fit_data(self,
            handle_limit=DEFAULT_HANDLE_LIMIT,
            upper_limit=DEFAULT_UPPER_LIMIT):
        f_handle, f_upper = np.linspace(*handle_limit), np.linspace(*upper_limit)
        MMS_global = curvefit.MMS_rip_maker(handle_limit, upper_limit)
        return  MMS_global(np.append(f_handle, f_upper), 
            *self.rip_params()
        )

    def test_return_equivalent_MMS_extension(self):
        global_fit_ext = self._generate_fit_data()
        self.assertItemsEqual(self.rip_test_ext, global_fit_ext, 
            msg='Global fit does not equal MMS equivalent\n{}\n{}'.format(
                global_fit_ext, self.rip_test_ext))

    def test_global_fit_equals_MMS_single_rip(self):
        handle_limit = (1,11)
        upper_limit = (10,20)
        ext, forces = self._generate_test_regions(
            [handle_limit, upper_limit],
            Lc=[10])
        handle_x, upper_x = ext
        handle_f, upper_f = forces
        fit = curvefit.fitWLCrip(handle_x, handle_f, upper_x, upper_f)
        print handle_x
        print fit(handle_f)
        print handle_x - fit(handle_f)
        print fit
        self.assertItemsAlmostEqual(handle_x, fit(handle_f))
        # self.assertItemsAlmostEqual(upper_x, fit(upper_f))
        print fit.parameters
        print self.rip_params()
        self.assertItemsAlmostEqual(fit.parameters.values(),
            self.rip_params(),
            delta=0.001)

    def assertItemsAlmostEqual(self, first, second, delta=0.1):
        map(lambda x,y: self.assertAlmostEqual(x, y, delta=delta), first, second)

def setUp():
    global X,F,fit_params
    X = range(700,1100,10)
    F = range(0,21,1)
    fit_params = (30, 1150, 0, 1e6)
  
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
  
