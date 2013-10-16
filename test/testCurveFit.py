import numpy as np
from smbanalyze import curvefit
import unittest
from mock import MagicMock, patch

X = None
Y = None
params = (27, 1046, 0, 1200)

# something like this to check for masks converting ok
# clf()
# for mask,converted in zip(mask_limits, converted_masks):
#     fplot.plot(p.trap[mask], hold=True, style='x', markersize=5)
#     fplot.plot(p.trap[fit.mask][converted], hold=True, markersize=3)

class MMS_rip_test(unittest.TestCase):

    DEFAULT_HANDLE_LIMIT = (1,10)
    DEFAULT_UPPER_LIMIT = (13,20)

    def setUp(self):
        self.handle_params = (30, 1150, 0, 1200)
        self.rip_test_ext, self.rip_test_force = self._generate_test_data()

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

    def test_fitwlc_masks(self):
        pass
        # fit = curvefit.fitWLC_masks(self.rip_test_ext, self.rip_test_force,
        #    self.rip_test_ext>850)

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
        print handle_f
        print fit(handle_f)
        print handle_x - fit(handle_f)
        print fit
        self.assertItemsAlmostEqual(handle_x, fit(handle_f))
        print upper_x
        print fit(upper_f)
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


class TestFitRegions(unittest.TestCase):
    def setUp(self):
        self.x = np.arange(800,1000,10)
        self.y = np.linspace(1, 15, len(self.x))

    @patch('smbanalyze.curvefit.MMS_rip_region_maker')
    def test_accepts_single_region(self, func):
        return
        region = self.x>900
        FitRegions2(self.x, self.x, [region])

    def test_accepts_two_regions(self):
        pass