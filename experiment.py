__all__ = ['split_reverse_pull', 'split_pulls_at_point', 
            'Pulling', 'ExpList', 'group_by']

import os.path as opath
from operator import itemgetter, attrgetter, methodcaller, and_
import operator as op
import logging
import cPickle as pickle
from re import search as re_search
from itertools import groupby, ifilter
from smbanalyze.date import today, date, to_date
from smbanalyze.sampling import downsample as dsample, down_average
from smbanalyze.pulling import Pulling
from smbanalyze.numlist import numlist

from matplotlib.mlab import find
from numpy import mean, all, any, diff, where, vstack

import fileIO 
import constants
import image
from image import ROI
from datatypes import TrapData, FretData, Data
from fancydict import nesteddict


logger = logging.getLogger(__name__)
logger.setLevel(constants.logLevel)
logger.addHandler(constants.logHandler)


def to_fret_ext(p, sampling_ratio=None):
  '''Returns (fret,ext) with downaveraged extension from experiment p
  sampling_ratio is read from p if it exists, but can be specified/overridden
  '''
  try:
    if sampling_ratio is None:
      sampling_ratio = p['sampling_ratio']
  except KeyError:
    raise Exception('sampling_ratio not found in experiment. Please specify!')
  return vstack( (p.fret.fret, down_average(p.trap.ext, sampling_ratio)) )

def loadexp(filename):
  basename, extension = opath.splitext(filename)
  filename = basename + (extension or '.exp')
  return pickle.load(open(filename,'rb'))

def fromMatch(*globs):
  flist = exp_files(*globs)
  cls = Pulling
  matched = filter(lambda fname: cls.filenameMatchesType(fname), flist)
  files_not_loading = set(flist)-set(matched)
  if Options.loading.filename_matching and files_not_loading:
    msg = "Files matching pattern were skipped because of option 'filename_matching' is ON:\n{}"
    logger.warning(msg.format('\n'.join(files_not_loading)))
    return ExpList(map(cls.fromFile, matched))
  else:
    return ExpList(map(cls.fromFile, flist))

def exp_files(*fglob):
  'Return list of unique filenames (without extension) of pulling and fret file types'
  return fileIO.files_matching(fglob, with_ext=(fileIO.PULL_FILE,),
    keep_ext=False)

def fromFiles(*filelist):
  'Load experiments from a list of files or an argument list'
  filelist = collapseArgList(filelist)
  return ExpList(map(fromFile, filelist))

def fromFile(filename):
  basename, ext = fileIO.splitext(filename)
  if not ext or Pulling.EXTENSION == ext:
    return Pulling.fromFile(filename)

def collapseArgList(arglist):
  'Allows a function to take either *args or a list as its argument'
  if isinstance(arglist[0], list) or isinstance(arglist[0], tuple):
    if len(arglist)>1: raise ValueError('Argument must be a list or *args')
    return arglist[0]
  else:
    return arglist

def group_by(iterable, keyfunc):
    return {key: ExpList(p) for key,p in groupby(iterable, keyfunc)}

def on_metadata(key):
  return lambda p: p[key]

def split_pulls_at_point(exps, point):
  '''Return tuple (below,above) distinguished by their relative position above/below "point"''' 
  ext_cutoff, force_cutoff = point
  high, low = ExpList(), ExpList()
  for p in exps:
    f_at_ext = p.trap.at(ext=ext_cutoff).f
    if not f_at_ext or f_at_ext > force_cutoff:
      high += [p]
    else:
      low += [p]
  return low, high

def loadimg(p, directory='.', **kwargs):
  '''Return the Stack image using ROI and background stored in experiment metadata.
  Must be in current or specified directory.
  '''
  filename = p.filename if directory=='.' else opath.join(directory, p.filename)
  meta = p.fret.metadata
  if 'roi' not in kwargs:
    # If keyword roi= is not specified,
    # get the ROIs from the experiment metadata
    rois = [ meta[roi_key] 
      for roi_key in ('roi_donor', 'roi_acceptor')
      if roi_key in meta
      ]
    kwargs['roi'] = map(ROI.fromDict, rois)
  kwargs.setdefault('background', meta.get('background', '') )
  return image.fromFile(fileIO.add_img_ext(filename), **kwargs)


class ExpList(numlist):
  def matching(self, *match):
    "Return ExpList with experiment names matching *match globs"
    return self.filter(lambda x: re_search(fileIO.makeMatchStrFromArgs(*match), x.filename))

  def has_meta(self, key, value=None):
    if value is not None:
      return self.filter(lambda p: p.metadata.get(key,None)==value)
    else:
      return self.filter(lambda p: key in p.metadata)
      
  def has_attr(self, *attributes):
    "Returns ExpList with elements that have given attributes != None"
    condition = lambda p: all(map(lambda name: getattr(p, name, None) is not None,
                          attributes))
    return self.filter(condition)

  def not_has(self, attr):
    "Returns ExpList of experiments which DON'T HAVE the given attribute set"
    condition = lambda p: getattr(p, attr, None) is None
    return self.filter(condition)

  def is_a(self, kind):
    "Returns ExpList of experiments of specified kind"
    return self.filter(lambda p: isinstance(p, kind))

  HAS_VALUE_COMPARE = {'greaterthan': op.gt, 'atleast': op.ge, 
                      'lessthan': op.lt, 'atmost': lambda x,y: not any(x>y),
                      'equals': op.eq}

  def has_value(self, **kwargs):
    '''Return ExpList with item matching filter arguments
    Options: greaterthan, atleast, lessthan, atmost, equals
    Examples:
    pulls.has_value(trap_f_atleast=15)
    pulls.has_value(trap_ext_greaterthan=750)
    pulls.has_value(fret_time_atleast=30)
    '''
    key,value = kwargs.popitem()
    def filter_by_value(to_filter, key, value):
      attr, comparison = key.rsplit('_', 1)
      attr = attr.replace('_','.')
      comparator = ExpList.HAS_VALUE_COMPARE.get(comparison, None)
      if comparator:
        return to_filter.filter( lambda x: any(comparator(attrgetter(attr)(x), value)) )
      else:
        raise ValueError('Comparison operator <{}> is not defined'.format(comparison))

    return filter_by_value(self, key, value)

  def collapse(self, trap_sorted_by='ext', fret_sorted_by='time'):
    "Collapse experiments into a single experiment with all data appended together."
    assert isinstance(trap_sorted_by, str)
    assert isinstance(fret_sorted_by, str)
    if not self._all_elements_have_attr('trap'):
      raise AttributeError(
        'All experiments in ExpList must have attribute "trap"')
    filtered_by_fret = self.has_attr('fret')
    num_with_fret = len(filtered_by_fret)
    fret_data = None
    if num_with_fret == len(self):
        fret_data = FretData.aggregate(self.getattrs('fret'), fret_sorted_by)
    elif num_with_fret > 0:
        logger.warning('Not all experiments have fret: not collapsing fret data!')
    trap_data = TrapData.aggregate(self.getattrs('trap'), sort_by=trap_sorted_by)
    fname = self[0].filename or ''
    fname += '_collapsed' if fname else 'collapsed'
    return Pulling(trap_data, fret_data, dict(filename=fname, collapsed=True))

  def _all_elements_have_attr(self, attr):
    for x in self:
      if not hasattr(x, attr):
        return False
    return True

  def adjustForceOffset(self, baseline=0.0, offset=None, x_range=None):
    return self.call('adjustForceOffset', baseline, offset, x_range)

  def adjustExtensionOffset(self, baseline=None, x_range=None, f_range=None):
    baseline = baseline or mean(self.call('extensionOffset',
                  x_range=x_range, f_range=f_range))
    return self.call('adjustExtensionOffset', baseline, 
                  x_range=x_range, f_range=f_range)

  def adjustOffset(self, to_f=0.5, to_x=None, 
          ext_f_range=None, ext_x_range=None,
          force_x_range=None):
    '''Adjust force and extension offsets returning xoffset, foffset
    
    to_f:  The force baseline that traces will be adjusted to. Default: 0.5 pN

    to_x:  The extension baseline. None (default) calculates average extension over
            f_range.
            
    ext_f_range: The range of forces used to calculate the extension offset.
            None (default) uses value Pulling.extensionOffsetRange

  ext_x_range: The range of extensions used to calculate the extension offset.
  
    force_x_range: The range of extensions used to calculate the force offset.
            None (default) uses value Pulling.forceOffsetRange
    '''
    if ext_f_range is not None:
        filter_f = ext_f_range[1]
        to_adjust = self.has_value(trap_f_atleast=filter_f)
    else:
        to_adjust = self.auto_filter()
        filter_f = Options.filtering.required_pulling_force
    if to_adjust != self:
      msg = 'Ignoring experiments below {} pN!'.format(filter_f)
      logger.warning(msg)
    foffset = to_adjust.adjustForceOffset(to_f, x_range=force_x_range)
    xoffset = to_adjust.adjustExtensionOffset(to_x, x_range=ext_x_range, f_range=ext_f_range)
    return xoffset, foffset

  def split_reverse_pull(self):
    return ExpList.map(split_reverse_pull, self.is_a(Pulling))

def split_reverse_pulls(exps):
  return ExpList.map(split_reverse_pull, exps.is_a(Pulling))

def split_reverse_pull(exp):
  '''
  fwd, rev = split_reverse_pull(exp)
  Return forward and reverse pulls as separate experiments
  '''
  assert isinstance(exp, Pulling)
  split = find_reverse_splitpoint(exp.trap)
  if split is None:
    return exp, None
  fwd, rev = split_reverse_data(exp.trap, split)
  fwd_fret, rev_fret = None, None
  if exp.fret:
    ratio = exp.metadata.get('sampling_ratio', 1)
    fwd_fret, rev_fret = split_reverse_data(exp.fret, split/ratio)
  rev_meta = exp.metadata.copy()
  rev_meta['reverse'] = True
  return (Pulling(fwd, fwd_fret, exp.metadata), 
          Pulling(rev, rev_fret, rev_meta)
        )

def split_reverse_data(data, split):
  '''Return forward and reverse data on split'''
  assert isinstance(data, Data)
  if split is None:
    raise ValueError('Split cannot be None')
  cls = type(data)
  forward, reverse = (cls.fromObject(data[:split]), 
          cls.fromObject(data[:split-1:-1]))
  reverse.metadata['reverse'] = True
  return forward, reverse

def find_reverse_splitpoint(trap):
  '''Return index location of reverse pull splitpoint or None if not found'''
  ext, force, sep = trap
  i = where(diff(sep)<=0)[0]
  if len(i) == 0:
    return None
  else:
    return i[0]



class OptionError(Exception):
  pass


class OptionDict(dict):
  def __init__(self, mapping):
    for key, value in mapping.items():
      if isinstance(value, dict):
        mapping[key] = OptionDict(value)
    super(OptionDict, self).__init__(mapping)

  def __getattr__(self, key):
    try:
      return self[key]
    except KeyError:
      raise OptionError('Option "{}" does not exist'.format(key))

  def __setattr__(self, key, val):
    if isinstance(val, dict) or not self.has_key(key):
      raise ValueError('Options are immutable! Can only set *existing* values.')
    self[key] = val

  def __repr__(self):
    out = []
    for key, value in self.items():
      if isinstance(value, dict):
        out.append("{}:\n\t{}".format(key, repr(value).replace('\n','\n\t')))
      else:
        out.append("{} = {}".format(key, value))
    return '\n'.join(out)


Options = OptionDict({
  'loading': {
    'filename_matching': False},
  'filtering': {
    'auto_filter': True, # apply these options automagically where needed
    'required_pulling_force': Pulling.extensionOffsetRange[1]}
  
})
