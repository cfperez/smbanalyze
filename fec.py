'''
'''

import numpy as np
from numpy import NAN
from collections import Iterable
from smbanalyze.curvefit import Fit,FitRegions
from useful import broadcast
from constants import kT, parameters

RIP_NAME_PREFIX = 'Lc'

def findrips(trap, min_rip_ext):
    handle_data = trap.select(x=(None,min_rip_ext))

    # the difference derivative of the force below min_rip_ext is
    # used as the baseline/expected derivative for WLC curve
    handle_deriv = np.diff(handle_data.f)
    
    # where the derivative first exceeds the minimum derivative found in
    # the handle region, call that the rip
    rip_location = np.find(np.diff(trap.f) < min(handle_deriv))
    if len(rip_location)==0:
        return np.asarray([NAN,NAN,NAN])
    rip_location = rip_location[0]
    return trap[rip_location]

# TODO: Currently redundant with TrapData.make_masks, which to use??
def make_masks(trap, intervals):
  return map(trap.mask_from_interval, intervals)

def rip_sizes(fit):
  assert isinstance(fit, Fit)
  lc = [fit[param]
    for param in fit if param.startswith(RIP_NAME_PREFIX)]
  # Replace first 'Lc' (handle contour) with 0 for
  # pairwise subtraction below
  lc[0] = 0
  return [lc[n]-lc[n-1] for n in range(1,len(lc))]

def rip_errors(fit):
  assert isinstance(fit, Fit)
  lc_err = [val for p,val in fit.error.items()
    if p.startswith(RIP_NAME_PREFIX) and p!=RIP_NAME_PREFIX]
  return lc_err[:1] + [np.sqrt(lc_err[n]**2+lc_err[n+1]**2)
    for n in range(len(lc_err)-1)]

DEFAULT_NA_TYPE = 'RNA'
RISE_PER_BASE = {'RNA': 0.59, 'DNA': 0.58}
HELIX_SIZE = {'RNA': 2.2, 'DNA': 2.0}
PERSISTENCE_LENGTH = dict(RNA=1.0, DNA=1.2)

def nm_to_nt(nm, stems_lost=1, helix_size=None, na_type='RNA'):
    '''Return basepairs from contour length change in nm assuming RNA
    stems_lost = # of helix stems lost during rip. <0 more helix stems exposed.
    '''
    helix_size = helix_size or HELIX_SIZE[na_type]
    return (nm+helix_size*stems_lost)/RISE_PER_BASE[na_type]

class Rips(object):
  DEFAULT_HELIX_SIZE = 2.2

  def __init__(self, rips, stems_lost, na_type, error=[]):
    assert isinstance(na_type, str) and na_type in HELIX_SIZE
    assert isinstance(rips, Iterable)
    if len(rips) != len(stems_lost):
      raise ValueError(
        "Number of stems_lost ({}) must equal number of rips ({})".format(
          stems_lost, rips))
    self._in_nm = rips
    self.na_type = na_type
    self.helix_size = HELIX_SIZE[na_type]
    self.stems_lost = stems_lost
    self._error = error

  @classmethod
  def from_fit(cls, fit, stems_lost, na_type):
    return cls(rip_sizes(fit), stems_lost, na_type, rip_errors(fit))

  @property
  def size_nm(self):
    return self._in_nm

  @property
  def size_nt(self):
    return map(
      lambda rip,stems_lost: nm_to_nt(rip, stems_lost, 
      helix_size=self.helix_size, na_type=self.na_type),
      self._in_nm,
      self.stems_lost
      )

  @property
  def error_nm(self):
    return self._error

  @property
  def error_nt(self):
    return map(
      lambda error: nm_to_nt(error, 0,
        helix_size=self.helix_size, na_type=self.na_type),
      self._error)

  def __repr__(self):
    return "Rips(rips={}, stems_lost={}, na_type='{}', error={})".format(
      self.size_nm, self.stems_lost, self.na_type, self._error)

  def __unicode__(self):
    if self._error:
      rip_size_fmt =  u'{:.1f}\u00B1{:.1f} nt ({:.2f} nm)'
      rips_str = "Rips: " + " | ".join(rip_size_fmt.format(nt,err,nm)
        for nt,err,nm in zip(self.size_nt,self.error_nt,self.size_nm))
      total_str = "Total: "+rip_size_fmt.format(
        sum(self.size_nt),self.error_nt[-1],sum(self.size_nm))
    else:
      rip_size_fmt = u'{:.1f} nt ({:.2f} nm)'
      rips_str = "Rips: " + " | ".join(rip_size_fmt.format(nt,nm)
        for nt,nm in zip(self.size_nt,self.size_nm))
      total_str = "Total: "+rip_size_fmt.format(sum(self.size_nt),sum(self.size_nm))
    return rips_str+'\n'+total_str

  def __str__(self):
    return unicode(self).encode('utf-8')

def analyze_rips(trap, intervals, stems_lost, 
  handle_above=None, na_type=None, **fitOptions):
  assert isinstance(stems_lost, (list, tuple))
  assert isinstance(intervals, (list, tuple))
  na_type = na_type or DEFAULT_NA_TYPE
  if na_type not in HELIX_SIZE:
    raise ValueError('na_type must be one of: {}'.format(HELIX_SIZE.keys()))
  masks = trap.make_masks(intervals)
  if handle_above:
    above = trap.mask_above(handle_above)
    masks[0] = masks[0] & above
  fitOptions.setdefault('Lp1', PERSISTENCE_LENGTH[na_type])
  fit = fit_rips(trap, masks, **fitOptions)
  fit.plot()
  return Rips.from_fit(fit, stems_lost, na_type), fit

def analyze(exp, intervals, stems_lost,
  handle_above=None, na_type=None, **fitOptions):
  rips,fit = analyze_rips(exp.trap, intervals,
       stems_lost, handle_above, na_type, **fitOptions)
  exp.rips = rips
  exp.fit = fit
  return rips

def fit_rips(trap, masks, **fitOptions):
  return fit_masks(trap.ext, trap.f, masks, **fitOptions)

def fit_masks(x, f, masks, **fitOptions):
  '''
  Return Fit using MMS_rip_region fit function on data in x and f.
  masks: list of boolean arrays which indicate:
    1) the handle (1st element)
    2) the rips (following elements)
  '''
  assert isinstance(x, Iterable) and len(x)>0
  assert isinstance(f, Iterable) and len(f)==len(x)
  assert map(lambda m: isinstance(m, np.ndarray), masks)
  assert map(lambda m: m.dtype is np.dtype('bool'), masks)

  fitOptions.setdefault('Lc', max(x[masks[-1]]))
  fitOptions.setdefault('fixed', tuple())
  fitOptions['fixed'] += ('K', 'K1', 'Lp1')
  fit = FitRegions(x, f, fitfunc=MMS_rip_region_maker,
    regions=masks,
    **fitOptions)
  return fit

def fit_wlc(trap, mask=None, fitfunc='MMS', **fitOptions):
  "Fit stretching data to WLC model with MMS default fitfunc"
  assert mask is None or len(mask) == len(trap)
  fitOptions.setdefault('Lc', max(trap.ext)*1.05)
  if isinstance(fitfunc, str):
    fitfunc = eval(fitfunc)
    fixed = getattr(fitfunc, 'fixed', None)
    if fixed:
      fitOptions.setdefault('fixed', fixed)
  fit = Fit(trap.ext, trap.f, mask=mask, fitfunc=fitfunc, **fitOptions)
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
MMS.fixed = 'K'

def moroz_nelson(F, Lp, Lc, F0, K):
  F_ = np.asarray(F)-F0
  return Lc * (1 - 0.5*pow(F_*Lp/kT(parameters['T'])-1/32.,-.5) + F_/K)
moroz_nelson.default = {'Lp':30., 'Lc':1150., 'F0':0.1, 'K': 1200.}
moroz_nelson.inverted = True
moroz_nelson.fixed = 'K'

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
    handle_ext = moroz_nelson(handle_force, Lp, Lc, F0, K)
    # Try to use moroz_nelson fitting function here...
    # Make into more modular form to allow comparisons
    # by passing in a fitfun= param?
    # rip_ext = [moroz_nelson(force,Lp,Lc,F0,K)+moroz_nelson(force,Lp1,Lc_rip,F0,K1)
    #             for force,Lc_rip in rip_items]
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
