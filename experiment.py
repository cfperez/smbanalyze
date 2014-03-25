__all__ = ['split_reverse_pull']

import os.path as opath
from operator import itemgetter, attrgetter, methodcaller
import operator as op
import logging
import cPickle as pickle
import re
import abc
from functools import total_ordering
from itertools import groupby
from datetime import datetime


from matplotlib.mlab import find
from numpy import min, max, asarray, insert, sum, mean, all, any, diff, NAN, where

import fileIO 
from curvefit import fitWLC, fitWLC_masks
from useful import groupat
import fplot
import constants
import image
from datatypes import TrapData, FretData, hasTrapData, hasFretData, AbstractData

logger = logging.getLogger(__name__)
logger.setLevel(constants.logLevel)
logger.addHandler(constants.logHandler)

PULL_FILENAME_INFO = fileIO.MOL_FILENAME_INFO + ('pull',)

###############

def aggregate(*data):
  return ''

###############


def mol_info_dict(*args):
  assert len(args) == len(PULL_FILENAME_INFO)
  return dict(zip(PULL_FILENAME_INFO, args))

class ExperimentError(Exception):
  pass

def load(filename):
  basename, extension = opath.splitext(filename)
  filename = basename + (extension or '.exp')
  return pickle.load(open(filename,'rb'))

def fromData(*datalist, **kwargs):
  "List of experiments from PullData and FretData type"
  exptype = kwargs.get('type','trap')
  if exptype=='trap':
    output = []
    for pull,fret in groupat(hasTrapData, datalist, size=2):
      output += [Pulling(pull,fret)]
    return output if len(output)>1 else output[-1]

def fromMatch(*globs):
  flist = exp_files(*globs)
  cls = Pulling
  matched = filter(lambda fname: cls.filenameMatchesType(fname), flist)
  files_not_loading = set(flist)-set(matched)
  if Options.loading.filename_matching and files_not_loading:
    msg = "Files matching pattern were skipped because of option 'filename_matching' is ON:\n{}"
    logger.warning(msg.format('\n'.join(files_not_loading)))
    return List(map(cls.fromFile, matched))
  else:
    return List(map(cls.fromFile, flist))

def exp_files(*fglob):
  'Return list of unique filenames (without extension) of pulling and fret file types'
  return fileIO.files_matching(fglob, with_ext=(fileIO.PULL_FILE, fileIO.FRET_FILE),
    keep_ext=True)

def fromFiles(*filelist):
  'Load experiments from a list of files or an argument list'
  filelist = collapseArgList(filelist)
  return List(map(fromFile, filelist))

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


class List(list):
  '''
  Ordered (and sorted) container for experiments
  '''
  def __init__(self, iterable=[]):
    super(List, self).__init__(iterable)
    self._figure = fplot.Figure()
    try:
      self.sort()
    except AttributeError:
      pass

  @property
  def figure(self):
    return self._figure

  def filter(self, condition):
    "Returns a List with experiments matching condition"
    return List( filter(condition, self) )

  def auto_filter(self):
    "Returns list with auto_filter options applied"
    options = Options.filtering
    autoforce = options.required_pulling_force
    return self.has_value(trap_f_atleast=autoforce) if autoforce else self

  def matching(self, *match):
    "Return List with experiment names matching *match globs"
    return self.filter(lambda x: re.search(fileIO.makeMatchStrFromArgs(*match), x.filename))

  def get(self, name, *more):
    "Get all attributes of experiments in List by name: mol.get('f')"
    try:
      return map(attrgetter(name, *more), self)
    except AttributeError:
      return map(itemgetter(name, *more), self)
    else:
      raise ExperimentError('Missing attribute {0} in a List element'.format(name))

  def has(self, *attributes):
    "Returns List of experiments which have all given attributes != None"
    condition = lambda p: all(map(lambda name: getattr(p, name, None) is not None,
                          attributes))
    return self.filter(condition)

  def not_has(self, attr):
    "Returns List of experiments which DON'T HAVE the given attribute set"
    condition = lambda p: getattr(p, attr, None) is None
    return self.filter(condition)

  def is_a(self, kind):
    "Returns List of experiments of specified kind"
    return self.filter(lambda p: isinstance(p, kind))
      
  def call(self, action, *args, **kwargs):
    "Call function <action> with *args and **kwargs on all experiments"
    try:
      return map( methodcaller(action, *args, **kwargs), self )
    except AttributeError:
      raise ExperimentError('Missing method {0} in a List element'.format(action))

  HAS_VALUE_COMPARE = {'greaterthan': op.gt, 'atleast': op.ge, 
                      'lessthan': op.lt, 'atmost': lambda x,y: not any(x>y),
                      'equals': op.eq}

  def has_value(self, **kwargs):
    '''
    See HAS_VALUE_COMPARE for list of operators
    Examples:
    pulls.has_value(trap_f_atleast=15)
    pulls.has_value(trap_ext_greaterthan=750)
    pulls.has_value(fret_time_atleast=30)
    '''
    for key, value in kwargs.items():
      attr, comparison = key.rsplit('_', 1)
      attr = attr.replace('_','.')
      comparator = List.HAS_VALUE_COMPARE.get(comparison, None)
      if comparator:
        return self.filter( lambda x: any(comparator(attrgetter(attr)(x), value)) )
      else:
        raise ValueError('Comparison operator <{}> is not defined'.format(comparison))

  def collapse(self, trap_sorted_by='ext', fret_sorted_by='time'):
    "Collapse experiments into a single experiment with all data appended together."
    assert isinstance(trap_sorted_by, str)
    assert isinstance(fret_sorted_by, str)
    if not self._all_elements_have_attr('trap'):
      raise ExperimentError(
        'All experiments in List must have attribute "trap"')
    filtered_by_fret = self.has('fret')
    num_with_fret = len(filtered_by_fret)
    fret_data = None
    if num_with_fret == len(self):
        fret_data = FretData.aggregate(self.get('fret'), fret_sorted_by)
    elif num_with_fret > 0:
        logger.warning('Not all experiments have fret: not collapsing fret data!')
    trap_data = TrapData.aggregate(self.get('trap'), sort_by=trap_sorted_by)
    fname = self[0].filename or ''
    fname += '_collapsed' if fname else 'collapsed'
    return Pulling(trap_data, fret_data, filename=fname, collapsed=True)

  def _all_elements_have_attr(self, attr):
    filtered_by_attr = self.has(attr)
    return len(filtered_by_attr) == len(self)

  def _aggregate_data(self, attr, sort_by):
    filtered_by_attr = self.has(attr)
    num_with_attr = len(filtered_by_attr)
    data = None
    if num_with_attr == len(self):
      data = AbstractData.aggregate(self.get(attr), sort_by)
    elif num_with_attr > 0:
      logger.warning(
        'Not all experiments have {atrr}: not collapsing {attr} data!'.format(attr=attr))
    return data

  def plot(self, **options):
    "Plot all experiments overlayed onto same panel. Can plot only trap or fret."
    options.setdefault('legend', None)
    self.figure.show().clear()
    for p in self:
      fplot.plot(p.fret, p.trap, 
        hold=True, 
        label=p.metadata.get('filename', ''),
        **options)
    return self.figure

  def savefig(self, filename):
    "Save figure from last plotall()"
    self.figure.toFile(filename)

  def saveallfig(self, path=''):
    "Save all individual plots"
    self.call('savefig', path=path)

  def fitHandles(self, *args, **kwargs):
    return self.call('fitHandles', *args, **kwargs)

  def findRip(self, min_rip_ext=None):
      return asarray(self.call('findRip', min_rip_ext))

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

  def save(self, filename):
    with open(filename, 'wb') as fh:
      pickle.dump(self, fh, protocol=2)

  def split_reverse_pull(self):
    return List(map(split_reverse_pull, self.is_a(Pulling)))

  def __getslice__(self, i, j):
    return self.__getitem__(slice(i, j))

  def __getitem__(self, key):
    out = super(List, self).__getitem__(key)
    if isinstance(out, list):
      return List(out)
    else:
      return out

  def next(self):
    try:
      return self.it.next()
    except (StopIteration, AttributeError):
      self.it = iter(self)
      return self.next()

  def __add__(self, other):
    return List(super(List, self).__add__(other))

  def __repr__(self):
    fmt = '\n '.join('{0}: {1}'.format(i, repr(item))
                      for i, item in enumerate(self))
    return '[' + fmt + ']'

def split_reverse_from_list(exps):
  return List(map(split_reverse_pull, exps))

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
  return (Pulling(fwd, fwd_fret, **exp.metadata), 
          Pulling(rev, rev_fret, reverse=True, **exp.metadata)
  )

def split_reverse_data(data, split):
  '''Return forward and reverse data on split'''
  assert isinstance(data, AbstractData)
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

def group_by(iterable, keyfunc):
    return {key: List(p) for key,p in groupby(iterable, keyfunc)}

def on_metadata(key):
  return lambda p: p.metadata.get(key, None)

@total_ordering
class Base(object):
  ".fret .f .ext and other meta-data (sample rate, trap speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  __metaclass__ = abc.ABCMeta

  INFO_FIELDS = ()

  def __init__(self, trap, fret, **metadata):
    if trap and not hasTrapData(trap):
      raise ExperimentError(
          "__init__ argument 'trap' <{}> does not have trap data".format(trap))
    if fret and not hasFretData(fret):
      raise ExperimentError(
          "__init__ argument 'fret' <{}> does not have fret data".format(fret))
    self.trap = trap
    self.fret = fret
    if trap:
      metadata.update(
        {('trap.' if not k.startswith('fret.') else '')+k: v 
        for k,v in trap.metadata.items() })
    if fret:
      metadata.update(
        {('fret.' if not k.startswith('trap.') else '')+k: v
        for k,v in fret.metadata.items() })
    self.metadata = metadata
    self._figure = fplot.Figure(self.filename)

  @abc.abstractmethod
  def filenameMatchesType(cls, filename):
    pass

  @property
  def filename(self):
    return self.metadata.get('filename', None)

  @property
  def info(self):
    return self.metadata.get('filename', '').split('_')
    return tuple(self.metadata.get(key, None) for key in self.INFO_FIELDS)

  def __lt__(self, other):
    this, that = self.info, other.info
    return this < that

  def __eq__(self, other):
    return self.info == other.info

  def __getitem__(self, key):
    return self.metadata[key]

  @property
  def figure(self):
    return self._figure

  def __setstate__(self, state):
    self.__dict__ = state
    self._figure = fplot.Figure(self.filename)

  @classmethod
  def fromFile(cls, strfile, fretfile, **metadata):
    assert strfile or fretfile
    assert isinstance(strfile, str) or strfile is None
    assert isinstance(fretfile, str) or fretfile is None

    trap = TrapData.fromFile(strfile) if strfile else None
    fret = FretData.fromFile(fretfile) if fretfile else None
    filename = strfile if strfile else fretfile
    if trap:
      today = datetime.today()
      trap.metadata.setdefault('date', (today.year-2000,today.month,today.day))
      metadata.update({'trap.'+k: v for k,v in trap.metadata.items()
        if not k.startswith('fret.')})
    metadata.setdefault('date', 
      trap.metadata.get('date', None) if trap else None)
    if fret:
      metadata.update({'fret.'+k: v for k,v in fret.metadata.items()
        if not k.startswith('trap.')})

    newCls = cls(trap, fret,
      filename=fileIO.splitext(filename)[0], 
      **metadata)
    assert isinstance(newCls, cls)
    assert hasattr(newCls, 'filename')
    assert newCls.filename is not None
    return newCls

  def save(self, filename=None, path=''):
    if filename is None and self.filename is None:
      raise ValueError('Please choose a filename')
    filename = filename or self.filename
    if not filename.endswith('.exp'):
      filename = filename + '.exp'
    full_path = opath.join(path, filename)
    with open(full_path, 'wb') as fh:
      pickle.dump(self, fh, protocol=2)

  @classmethod
  def fromMatch(cls, *filename_pattern):
    flist = filelist(*filename_pattern)
    matched = filter(lambda fname: cls.filenameMatchesType(fname), flist)
    files_not_loading = set(flist)-set(matched)
    if Options.loading.filename_matching and files_not_loading:
      msg = "Files matching pattern were skipped because of option 'filename_matching' is ON:\n{}"
      logger.warning(msg.format('\n'.join(files_not_loading)))
      return List(map(cls.fromFile, matched))
    else:
      return List(map(cls.fromFile, flist))

  def __repr__(self):
    if getattr(self, 'filename', None):
      return "<Experiment.%s from '%s'>" % (self.__class__.__name__, self.filename)
    else:
      return super(Base,self).__repr__()

  def __str__(self):
    return self.__repr__()

  @abc.abstractmethod
  def plot(self):
    pass

class Pulling(Base):
  "stretching curves of single molecule .str camera .cam image .img"

  forceOffsetRange = (740,770)
  extensionOffsetRange = (13,16)
  maximum_fitting_force = 25
  
  EXTENSION = fileIO.PULL_FILE
  INFO_FIELDS = ('date', 'construct', 'conditions', 'slide', 'mol', 'pull', 'reverse', 'collapsed' )

  def __init__(self, trap, fret=None, **metadata):
    super(Pulling, self).__init__(trap, fret, **metadata)
    self.handles = None
    self.lastFit = None
    self.fit = None
    self._ext_offset = 0

    if 'trap.sampling_time' in self.metadata and 'trap.step_size' in self.metadata:
      t,s = self.metadata['trap.sampling_time'], self.metadata['trap.step_size']
      self.metadata.setdefault('trap.pulling_rate', s/t)
    if fret:
      trap_rate = trap.metadata.get('sampling_time',
        constants.default_pulling_sampling_time) * 1000
      fret_rate = fret.metadata.get('exposurems',
        constants.default_fret_exposure_time_ms)
      if fret_rate % trap_rate != 0:
        raise ExperimentError(
          'Trap and FRET collection rates are not even multiples! {} (trap) vs {} (fret)'.format(
              trap_rate, fret_rate))
      self.metadata['sampling_ratio'] = int(fret_rate / trap_rate)

  @classmethod
  def parse_filename(self, basename):
    info, the_rest = fileIO.parse_mol_info(basename)
    pull = 1
    if the_rest:
      try:
          pull = int(the_rest)
      except:
        raise ValueError("Basename '%s' has incorrect format" % basename)
    info.update(pull=pull)
    return info

  @classmethod
  def fromFile(cls, strfile, fretfile=None, **metadata):
    'Load stretching data and corresponding fret data from files'
    basename, ext = fileIO.splitext(strfile)
    strfile, fretfileFromBase = (fileIO.add_pull_ext(basename),
                                fileIO.add_fret_ext(basename))
    fretfile = fretfile or fretfileFromBase
    if not opath.exists(fretfile):
      fretfile = None
    return super(Pulling, cls).fromFile(strfile, fretfile, **metadata)

  FILENAME_SYNTAX = ('construct', 'conditions', 'slide', 'mol', 'pull')

  @classmethod
  def filenameMatchesType(cls, filename):
    assert isinstance(filename, str)
    finfo = fileIO.parseFilename(filename)
    if not finfo:
      return False
    for attr in finfo._fields:
      val = getattr(finfo, attr)
      if attr not in cls.FILENAME_SYNTAX and val:
        return False
    return True

  @property
  def isReverse(self):
    '''Return True iff experiment is a reverse pull'''
    return self.metadata.get('reverse', False)

  def loadimg(self, directory='.', **kwargs):
    filename = self.filename if directory=='.' else opath.join(directory, self.filename)
    try:
      return image.fromFile(fileIO.add_img_ext(filename), **kwargs)
    except IOError:
      raise ExperimentError('IOError loading file: check image file location!')

  def forceOffset(self, x_range=None):
    x_range = x_range or Pulling.forceOffsetRange
    data = self.trap.select(x=x_range)
    if len(data) == 0:
      raise ExperimentError(
        '{0}: No data exists in range {1} - {2}'.format(
          str(self), *x_range))
    return mean(data.f)

  def adjustForceOffset(self, baseline=0.0, offset=None, offset_range=None):
    offset = offset or -self.forceOffset(offset_range)
    offset += baseline
    self.trap.f += offset
    self.trap.ext -= offset/self.trap.meanStiffness()
    return offset

  def extensionOffset(self, f_range=None, x_range=None):
    'Returns average extension of FEC between given forces'
    f_range = f_range or Pulling.extensionOffsetRange
    data = self.trap.select(x=x_range, f=f_range)
    if len(data) == 0:
      raise ExperimentError(
        '{0}: No data exists in f_range={1} and x_range={2}'.format(
          self, f_range, x_range)
      )
    return mean(data.ext)

  def adjustExtensionOffset(self, baseline, f_range=None, x_range=None):
    'Adjust extension to hit baseline. If f_range is not given, it is calculated from Pulling.extensionOffsetRange'
    offset = self.extensionOffset(x_range=x_range, f_range=f_range) - baseline
    self.trap.ext -= offset
    self._ext_offset = offset
    return offset

  # TODO: Move to TrapData, change name, and add docstring
  def recalculate(self, stiffness=None):
    if stiffness and len(stiffness) != 2:
      raise ValueError('Stiffness must be 2-tuple')
    current_k = self.trap.metadata.get('stiffness', constants.stiffness)
    new_k = stiffness or current_k
    self.trap.metadata['stiffness'] = new_k
    beadRadii = self.trap.metadata.get('bead_radii', constants.sumOfBeadRadii)

    displacement = self.trap.f/min(current_k)
    ratio = 1+min(new_k)/max(new_k)

    self.trap.f = displacement*min(new_k)
    self.trap.ext = self.trap.sep - beadRadii - displacement*ratio - self._ext_offset
    return self

  def plot(self, **kwargs):
    FEC = kwargs.setdefault('FEC', hasattr(self, 'fit') or not self.fret)
    if FEC:
      kwargs.setdefault('style', {'trap': '-'})
    #kwargs.setdefault('style', '-' if FEC else '.')
    kwargs.setdefault('legend', None)
    kwargs.setdefault('title', self.filename or '')
    kwargs.setdefault('label', self.filename or '')
    loc_x = min(self.trap.ext)+10
    location = list(kwargs.pop('annotate', (loc_x, 15)))
    if self.fret:
      self.figure.plot(self.fret, self.trap, **kwargs)
    else:
      self.figure.plot(self.trap, **kwargs)
    # tk if reverse pull, reverse the x-axis of the FEC, which is last thing plotted
    self.figure.xlim(reverse=self.isReverse)
    return self.figure

  def pickPoints(self, num=1):
    return self.figure.pickPoints(num)

  def pickRegions(self, num=1):
    return self.figure.pickRegions(num)
	
  def savefig(self, filename=None, path=''):
    if not self.figure.exists:
      raise ExperimentError('No figure available for experiment {0}'.format(str(self)))
    else: 
      filename = filename or self.filename
      self.figure.toFile(filename)

  @classmethod
  def load(cls, filename):
    basename, extension = opath.splitext(filename)
    filename = basename + (extension or '.exp')
    return pickle.load(open(filename,'rb'))


class OpenLoop(Base):
  "Object for manipulating FretData (and optional TrapData) for open loop measurements"
  def __init__(self, trap, fret, **metadata):
    super(OpenLoop, self).__init__(trap, fret, **metadata)

  FILENAME_SYNTAX = ('min', 'sec', 'force')

  @classmethod
  def filenameMatchesType(cls, filename):
    assert isinstance(filename, str)
    finfo = fileIO.parseFilename(filename)
    if not finfo:
      return False
    for attr in OpenLoop.FILENAME_SYNTAX:
      if getattr(finfo, attr) is not None:
        return True
    return False

  @classmethod
  def fromFile(cls, fretfile):
    assert isinstance(fretfile, str)
    basename, strfile, fretfile = fileIO.filesFromName(fretfile)
    if not opath.exists(strfile):
      strfile = None
    return super(OpenLoop, cls).fromFile(strfile, fretfile)

  def plot(self, **kwargs):
    self.figure.plot(self.fret, **kwargs)


class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass


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
