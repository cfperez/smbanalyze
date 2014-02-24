'''
'''

import numpy as np
from collections import Iterable
from smbanalyze.curvefit import FitRegions
from useful import broadcast
from constants import kT, parameters

RIP_NAME_PREFIX = 'Lc'

def make_masks(trap, intervals):
  return map(trap.mask_from_interval, intervals)

def rip_sizes(fit):
  lc = [fit[param] 
    for param in fit if param.startswith(RIP_NAME_PREFIX)]
  # Replace first 'Lc' (handle contour) with 0 for
  # pairwise subtraction below
  lc[0] = 0
  return [lc[n]-lc[n-1] for n in range(1,len(lc))]

DEFAULT_NA_TYPE = 'RNA'
RISE_PER_BASE = {'RNA': 0.59, 'DNA': 0.34}
HELIX_SIZE = {'RNA': 2.2, 'DNA': 2.0}

def nm_to_nt(nm, stems_lost=1, helix_size=2.2, na_type='RNA'):
    '''Return basepairs from contour length change in nm assuming RNA
    stems_lost = # of helix stems lost during rip. <0 more helix stems exposed.
    '''
    return (nm+helix_size*stems_lost)/RISE_PER_BASE[na_type]

class Rips(object):
  DEFAULT_HELIX_SIZE = 2.2

  def __init__(self, rips, stems_lost, na_type):
    assert isinstance(na_type, str) and na_type in HELIX_SIZE
    assert isinstance(rips, Iterable)
    self._in_nm = rips
    self.na_type = na_type
    self.helix_size = HELIX_SIZE[na_type]
    self.stems_lost = stems_lost

  @classmethod
  def from_fit(cls, fit, stems_lost, na_type):
    return cls(rip_sizes(fit), stems_lost, na_type)

  @property
  def in_nm(self):
    return self._in_nm

  @property
  def in_nt(self):
    return map(
      lambda rip,stems_lost: nm_to_nt(rip, stems_lost, 
      helix_size=self.helix_size, na_type=self.na_type),
      self._in_nm,
      self.stems_lost
      )

def analyze_rips(trap, intervals, stems_lost, 
  handle_above=None, na_type=None, **fitOptions):
  na_type = na_type or DEFAULT_NA_TYPE
  if na_type not in HELIX_SIZE:
    raise ValueError('na_type must be one of: {}'.format(HELIX_SIZE.keys()))
  masks = trap.make_masks(intervals)
  above = trap.mask_above(handle_above)
  masks[0] = masks[0] & above
  fit = fit_rips(trap, masks, **fitOptions)
  return Rips.from_fit(fit, stems_lost, na_type), fit

def fit_rips(trap, masks, **fitOptions):
  return fit_masks(trap.ext, trap.f, masks, **fitOptions)

def fit_masks(x, f, masks, **fitOptions):
  '''
  Return Fit using MMS_rip_region fit function on data in x and f.
  masks: list of boolean arrays which indicate:
    1) the handle (1st element)
    2) the rips (following elements)
  '''
  assert isinstance(x, Iterable)
  assert isinstance(f, Iterable)
  assert map(lambda m: isinstance(m, np.ndarray), masks)
  assert map(lambda m: m.dtype is np.dtype('bool'), masks)

  fitOptions.setdefault('Lc', max(x[masks[-1]]))
  fitOptions.setdefault('fixed', tuple())
  fitOptions['fixed'] += ('K', 'K1', 'Lp1')
  fit = FitRegions(x, f, fitfunc=MMS_rip_region_maker,
    regions=masks,
    **fitOptions)
  return fit

@broadcast
def MMS(F, Lp, Lc, F0, K):
  "Modified Marko-Siggia model as a function of force"
  f = float(F-F0)* Lp / kT(parameters['T'])
  inverted_roots = np.roots([1.0, f-0.75, 0.0, -0.25])
  root_index = int(f>=0.75)*2
  root_of_inverted_MS = np.real(inverted_roots[root_index])
  return Lc * (1 - root_of_inverted_MS + (F-F0)/float(K))
MMS.default = {'Lp':30., 'Lc':1150., 'F0':0.1, 'K': 1200.}
MMS.inverted = True

def MMS_rip(F, Lp, Lc, F0, K, Lp1, Lc1, K1):
  return MMS(F, Lp, Lc, F0, K) + MMS(F, Lp1, Lc1, F0, K1)
MMS_rip.default = dict(MMS.default, Lp1=1.0, Lc1=0., K1=1600.)
MMS_rip.inverted = True

def MMS_rip_region_maker(masks):
  '''Creates a fitting function for an arbitrary number of fit regions from masks
  that allows simultaneous fitting of each as MMS rips.
  '''
  assert isinstance(masks, (tuple,list))
  assert len(masks) > 1
  assert map(lambda e: isinstance(e, np.ndarray), masks)
  assert map(lambda e: e.dtype == np.dtype('bool'), masks)

  handle = masks.pop(0)

  def MMS_rip_region(F, Lp, Lc, F0, K, Lp1, Lc1, K1, **rips):
    rips['Lc1'] = Lc1
    handle_force = F[handle]
    rip_forces = [F[mask] for mask in masks]
    rip_sizes = map(lambda r: rips[r], sorted(rips))
    rip_items = zip(rip_forces, rip_sizes)
    handle_ext = MMS(handle_force, Lp, Lc, F0, K)
    rip_ext = [MMS_rip(force, Lp, Lc, F0, K, Lp1, Lc_rip, K1)
                for force,Lc_rip in rip_items]
    if len(rip_ext) > 1:
      rip_ext = np.concatenate(rip_ext)
    return np.append(handle_ext, rip_ext)
  
  addl_rips = ['Lc{}'.format(n) for n in range(2,1+len(masks))]
  MMS_rip_region.default = MMS_rip.default.copy()
  MMS_rip_region.default.update([(rip,10) for rip in addl_rips])
  MMS_rip_region.arglist = ['F', 'Lp', 'Lc', 'F0', 'K', 'Lp1', 'Lc1', 'K1'] + addl_rips
  MMS_rip_region.inverted = True
  lowest = ['Lp', 'Lc', 'F0', 'K']
  MMS_rip_region.args_by_region = [lowest] + [lowest+['Lp1', Lc, 'K1'] for Lc in ['Lc1']+addl_rips]

  return MMS_rip_region

#########
## Attempting to make MMS faster by using array
## Doesn't evaluate much faster than MMS using broadcast decorator
#########
def MMS2(F, Lp, Lc, F0, K):
  F_F0 = array(F)-F0
  f = F_F0* Lp / constants.kT(constants.parameters['T'])
  root_of_inverted_MS = _MMS_roots(f)
  return Lc * (1 - root_of_inverted_MS + F_F0/float(K))

def _root_index(f):
    return 2 if f>=.75 else 0

def _MMS_roots(f):
  return  array([ real(roots([1.0, f_-0.75, 0.0, -0.25])[_root_index(f_)]) for f_ in f])
