import os.path as path
from operator import itemgetter, attrgetter, methodcaller, add
import operator as op
import logging
import cPickle as pickle
import re
import abc

from matplotlib.mlab import find
from numpy import min, max, asarray, insert, sum, mean, all, any, diff, std, vstack, NAN
import matplotlib.pyplot as plt

import fileIO 
from curvefit import fitWLC
from useful import groupat
import fplot
import constants
import image
from datatypes import TrapData, FretData, hasTrapData, hasFretData

logger = logging.getLogger(__name__)
logger.setLevel(constants.logLevel)
logger.addHandler(constants.logHandler)

class ExperimentError(Exception):
  pass

def fromData(*datalist, **kwargs):
  "List of experiments from PullData and FretData type"
  exptype = kwargs.get('type','trap')
  if exptype=='trap':
    output = []
    for pull,fret in groupat(hasTrapData, datalist, size=2):
      output += [Pulling(pull,fret)]
    return output if len(output)>1 else output[-1]

def fromMatch(*globs):
  return Pulling.fromMatch(*globs)

def fromMatchAll(*fglob):
  'Load Pulling experiments with filenames that contain all of the strings listed in the arguments'
  fglob = list(fglob)
  files = filelist(*fglob)
  fglob[0] = path.basename(fglob[0])
  if not files:
    raise ExperimentError("No files found matching glob '{0}'".format(fglob))
  return fromFiles(files)
  
def filelist(*fglob):
  'Return list of unique filenames (without extension) of pulling and fret file types'
  files = fileIO.filtered_flist(*fglob, extensions=(fileIO.PULL_FILE, fileIO.FRET_FILE))
  files = map(lambda s: fileIO.splitext(s)[0], files)
  return list(set(files))

def fromFiles(*filelist):
  'Load experiments from a list of files or an argument list'
  filelist = collapseArgList(filelist)
  return List(map(fromFile, filelist))

def fromFile(filename):
  if OpenLoop.filenameMatchesType(filename):
    return OpenLoop.fromFile(filename)
  else:
    return Pulling.fromFile(filename)

def collapseArgList(arglist):
  'Allows a function to take either *args or a list as its argument'
  if isinstance(arglist[0], list) or isinstance(arglist[0], tuple):
    if len(arglist)>1: raise ValueError('Argument must be a list or *args')
    return arglist[0]
  else:
    return arglist

class List(list):
  def __init__(self, iterable=[]):
    super(List, self).__init__(iterable)
    try:
      self.sort(key=attrgetter('info.pull'))
      self.sort(key=attrgetter('info.mol'))
    except AttributeError:
      pass

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
      return map( attrgetter(name, *more), self)
    except AttributeError:
      raise ExperimentError('Missing attribute {0} in a List element'.format(name))

  def has(self, *attributes):
    "Returns List of experiments which have all given attributes not None"
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

  HAS_VALUE_COMPARE = {'greaterthan': op.gt, 'atleast': op.ge, 'lessthan': op.lt, 'atmost': op.le}

  def has_value(self, **kwargs):
    '''
    See HAS_VALUE_COMPARE for list of operators
    Examples:
    >>> pulls.has_value(trap_f_atleast=15)
    >>> pulls.has_value(trap_ext_greaterthan=750)
    >>> pulls.has_value(fret_time_atleast=30)
    '''
    for key, value in kwargs.items():
      attr, comparison = key.rsplit('_', 1)
      attr = attr.replace('_','.')
      try:
        comparator = List.HAS_VALUE_COMPARE[comparison]
        return self.filter( lambda x: any(comparator(attrgetter(attr)(x), value)) )
      except AttributeError:
        raise ValueError('Comparison operator <{}> is not defined'.format(comparison))

  METADATA_CHECK = {'trap': ('step_size', 'sampling_time'),
                    'fret': ('exposurems', 'frames', 'gain', 'binning') }

  def aggregate(self, datatype='trap', check_metadata=None):
    "Aggregate data of given kind from experiments by appending."
    if check_metadata is not None:
      assert isinstance(check_metadata, tuple) or isinstance(check_metadata, list)

    filtered = self.has(datatype)
    logger.info('Using experiments {}'.format(filtered))
    aggregated = reduce(add, filtered.get(datatype)) if filtered else None

    check_for_fields = check_metadata or List.METADATA_CHECK[datatype]
    for field in check_for_fields:
      try:
        meta = filtered.get(datatype+'.metadata')
        values = frozenset(map(itemgetter(field), meta))
        if len(values) > 1:
          logger.warning(
              'Multiple values found for metadata field <{}>: {}'.format(
               field, tuple(values))
               )
      except KeyError:
        logger.warning(
            'Not all items in List have metadata field <{}>'.format(
              field)
            )
    return aggregated

  def collapse(self, *types):
    "Collapse experiments into a single experiment with all data appended together."
    assert set(types).issubset(List.METADATA_CHECK.keys())
    fec = self.aggregate('trap')
    if fec is None:
      raise ExperimentError('No experiments in List had any trap data!')
    fdata = self.aggregate('fret')
    return Pulling(fec, fdata)

  def plotall(self, attr=None, **options):
    "Plot all experiments overlayed onto same panel. Can plot only trap or fret."
    assert attr in set([None,'trap','fret'])
    options.setdefault('labels', self.get('filename'))
    if hasattr(self,'figure') and self.figure.exists:
     self.figure.makeCurrent()
     self.figure.clear()
    if attr is None and self.has('fret'):
      options.setdefault('FEC', False)
      options.setdefault('legend', None)
      fplot.plotall(self.get('fret'), self.get('trap'), hold=True, **options)
    else:
      attr = attr or 'trap'
      fplot.plotall( self.get(attr), hold=True, **options)
    self.figure = fplot.Figure.fromCurrent()
    return self.figure

  def savefig(self, filename):
    "Save figure from last plotall()"
    self.figure.toFile(filename)

  def plot(self, **options):
    "Create individual plots"
    self.call('plot', **options)

  def saveallfig(self, path='.'):
    "Save all individual plots"
    self.call('savefig', path=path)

  def fitHandles(self, *args, **kwargs):
    return self.call('fitHandles', *args, **kwargs)

  def findRip(self, min_rip_ext=None):
      return asarray(self.call('findRip', min_rip_ext))
  
  def fitRip(self, *args, **kwargs):
    return self.call('fitRip', *args, **kwargs)

  def adjustForceOffset(self, baseline=0.0, offset=None, offset_range=None):
    return self.call('adjustForceOffset', baseline, offset, offset_range)

  def adjustExtensionOffset(self, baseline=None, offset_range=None):
    baseline = baseline or mean(self.call('extensionOffset', offset_range))
    return self.call('adjustExtensionOffset', baseline)

  def adjustOffset(self, to_f=0.5, to_x=None, f_range=None, x_range=None):
    '''Adjust force and extension offsets returning xoffset, foffset
    
    to_f  The force baseline that traces will be adjusted to. Default: 0.5 pN

    to_x  The extension baseline. None (default) calculates average extension over
            f_range.
            
    f_range The range of forces used to calculate the extension offset.
            None (default) uses value Pulling.extensionOffsetRange

    x_range The range of extensions used to calculate the force offset.
            None (default) uses value Pulling.forceOffsetRange
    '''
    if f_range is not None:
        filter_f = f_range[1]
        to_adjust = self.has_value(trap_f_atleast=filter_f)
    else:
        to_adjust = self.auto_filter()
        filter_f = Options.filtering.required_pulling_force
    if to_adjust != self:
      msg = 'Ignoring experiments below {} pN!'.format(filter_f)
      logger.warning(msg)
      self = to_adjust
    foffset = self.adjustForceOffset(to_f, offset_range=x_range)
    xoffset = self.adjustExtensionOffset(to_x, offset_range=f_range)
    return xoffset, foffset

  def saveall(self):
    self.call('save')

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

class RipAnalysis(object):
  '''Analyze rips
  rips = RipAnalysis.fromExperiments( List )
  RipAnalysis.fromTrapData( list ) ??
  rips.f
  rips.ext
  rips.sep
  mean rips.f
  std rips.f
  '''
  DEFAULT_MIN_RIP_EXT = 970
  
  def __init__(self, *exps):
    self.experiments = List()
    self.addExperiments(exps)
    
  @classmethod
  def fromExperiments(cls, *exps):
    newRips = cls(*exps)
    newRips.calculate(RipAnalysis.DEFAULT_MIN_RIP_EXT)
    return newRips

  def addExperiments(self, *exps):
    exps = collapseArgList(exps)
    self.experiments += List(exps)
    self.calculate(RipAnalysis.DEFAULT_MIN_RIP_EXT)

  def calculate(self, min_rip_ext=None):
    self.rips = self.experiments.findRip(min_rip_ext)
    
  def mean(self):
    return mean(self.rips, axis=0)

  def std(self):
    return std(self.rips, axis=0)

  def statistics(self):
    return vstack((self.mean(), self.std()))
    
  @property
  def f(self):
    return self.rips[:,1]

  @property
  def ext(self):
    return self.rips[:,0]
  
  def plot(self):
    raise NotImplementedError
  

class TrapSettings(object):
  pass

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, trap speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  __metaclass__ = abc.ABCMeta

  def __init__(self, trap, fret, **metadata):
    if trap and not hasTrapData(trap):
      raise ExperimentError(
          "__init__ argument 'trap' <{}> does not have trap data".format(trap))
    if fret and not hasFretData(fret):
      raise ExperimentError(
          "__init__ argument 'fret' <{}> does not have fret data".format(fret))
    self.trap = trap
    self.fret = fret
    self.metadata = metadata
    self.figure = fplot.Figure()
    self.filename = ''

  @abc.abstractmethod
  def filenameMatchesType(cls, filename):
    pass

  @classmethod
  def fromFile(cls, strfile, fretfile):
    assert strfile or fretfile
    assert isinstance(strfile, str) or strfile is None
    assert isinstance(fretfile, str) or fretfile is None

    trap = TrapData.fromFile(strfile) if strfile else None
    fret = FretData.fromFile(fretfile) if fretfile else None
    newCls = cls(trap, fret)
    filename = strfile if strfile else fretfile
    newCls.filename = fileIO.splitext(filename)[0]
    newCls.info = fileIO.parseFilename(newCls.filename)
    if not newCls.info:
      logger.warning('Problem parsing filename %s' % newCls.filename)
    assert isinstance(newCls, cls)
    assert hasattr(newCls, 'filename')
    assert newCls.filename is not None
    return newCls

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
    if getattr(self,'filename', None):
      return "<Experiment.%s from '%s'>" % (self.__class__.__name__, self.filename)
    else:
      return super(Base,self).__repr__()

  def __str__(self):
    return self.__repr__()

  def plot(self):
    raise NotImplementedError
    
class Pulling(Base):
  "stretching curves of single molecule .str camera .cam image .img"

  forceOffsetRange = (740,770)
  extensionOffsetRange = (13,16)

  def __init__(self, pull, fret=None, **metadata):
    super(Pulling, self).__init__(pull, fret, **metadata)
    self.handles = None
    self.lastFit = None
    self.resetRips()
    self._ext_offset = 0


  @classmethod
  def fromFile(cls, strfile, fretfile=None):
    'Load stretching data and corresponding fret data from files'
    basename, strfile, fretfileFromBase = fileIO.filesFromName(strfile)
    fretfile = fretfile or fretfileFromBase
    if not path.exists(fretfile):
      fretfile = None
    return super(Pulling, cls).fromFile(strfile, fretfile)

  FILENAME_SYNTAX = ('construct', 'conditions', 'slide', 'mol', 'pull')

  @classmethod
  def filenameMatchesType(cls, filename):
    assert isinstance(filename, str)
    finfo = fileIO.parseFilename(filename)
    if not finfo:
      logger.warning('Problem parsing filename "{}"'.format(filename))
      return False
    for attr in finfo._fields:
      val = getattr(finfo, attr)
      if attr not in cls.FILENAME_SYNTAX and val:
        return False
    return True

  def loadimg(self, directory='.', **kwargs):
    filename = self.filename if directory=='.' else path.join(directory, self.filename)
    try:
      return image.fromFile(fileIO.add_img_ext(filename), **kwargs)
    except IOError:
      raise ExperimentError('IOError loading file: check image file location!')

  def findRip(self, min_rip_ext=None):
    assert self.trap is not None
    min_rip_ext = self.handles.ext_range[1] if not min_rip_ext and self.handles \
        else min_rip_ext or None
    if min_rip_ext is None:
      raise ValueError('Must specify a min_rip_ext below which no rips occur')

    handle_data = self.trap.select(x=(None,min_rip_ext))

    # the difference derivative of the force below min_rip_ext is
    # used as the baseline/expected derivative for WLC curve
    handle_deriv = diff(handle_data.f)
    
    # where the derivative first exceeds the minimum derivative found in
    # the handle region, call that the rip
    rip_location = find(diff(self.trap.f) < min(handle_deriv))
    if len(rip_location)==0:
        return asarray([NAN,NAN,NAN])
    rip_location = rip_location[0]
    return self.trap[rip_location].data

  def fitHandles(self, x=None, f=None, **fitOptions):
    'Fit a WLC model to the lower part of the FEC corresponding to the handle stretching'
    fitOptions.setdefault('fitfunc','MMS')
    self.handles = self.fitForceExtension(x, f, **fitOptions)
    return self.handles

  def fitRip(self, x=None, f=None, guess=10, **fitOptions):
    'Fit a WLC model to the upper part of the FEC after the rip'
    if not self.handles:
      raise ExperimentError("Must fitHandles() before fitting rip!")
    parameters = self.handles.parameters.copy()
    fitOptions.setdefault('fitfunc', 'MMS_rip')
    parameters.update(fitOptions)
    parameters.setdefault('Lc1', guess)
    MMS_fixed = ('K','Lc','Lp','F0')
    if parameters.has_key('fixed'):
      parameters['fixed'] += MMS_fixed
    else:
      parameters['fixed'] = MMS_fixed + ('K1','Lp1')
    rip = self.fitForceExtension(x, f, **parameters)
    self.addRip(rip)
    return rip
    
  def addRip(self, newFit):
    if self._appendNewRip:
      self.rips += [newFit]
      self._appendNewRip = False
    else:
      self.rips[-1] = newFit

  def nextRip(self):
    self._appendNewRip = True

  def resetRips(self):
    self.rips = []
    self._appendNewRip = True

  @property
  def ripSizes(self):
    rips = asarray(map(itemgetter('Lc1'), self.rips))
    if len(rips) == 1:
      return rips
    else:
      rips = insert(rips, 0, 0)
      return diff(rips)

  def fitForceExtension(self, x=None, f=None, start=0, stop=-1, **fitOptions):
    "Fit WLC to Pulling curve and plot"
    pull = self.trap.select(x, f, (start,stop))
    if len(pull)==0:
      raise ValueError(
        'No trap data in interval defined by arguments x={} and f={}'.format(
          x,f))
    fit = fitWLC(pull.ext, pull.f, **fitOptions)
    fit.ext_range = (min(pull.ext), max(pull.ext))
    fit.f_range = (min(pull.f), max(pull.f))
    self.lastFit = fit
    if self.figure.exists:
      self.figure.plot(fit, hold=True)
    return fit

  @property
  def fits(self):
    if self.handles:
      return [self.handles]+self.rips
    else:
      return [self.lastFit] if self.lastFit else []

  def forceOffset(self, x_range=None):
    x_range = x_range or Pulling.forceOffsetRange
    data = self.trap.select(x=x_range)
    if len(data) == 0:
      raise ExperimentError(
        '{0}: No data exists in range {1} - {2}'.format(
          str(self), *x_range))
    return mean(data.f)

  def meanStiffness(self, stiffness=None):
    inverseAverage = lambda args: 1/sum(map(lambda x: 1./x, args))
    return inverseAverage(self.trap.metadata.get('stiffness', constants.stiffness))

  def adjustForceOffset(self, baseline=0.0, offset=None, offset_range=None):
    offset = offset or -self.forceOffset(offset_range)
    offset += baseline
    self.trap.f += offset
    self.trap.ext -= offset/self.meanStiffness()
    return offset

  def extensionOffset(self, frange=None):
    'Returns average extension of FEC between given forces'
    frange = frange or Pulling.extensionOffsetRange
    data = self.trap.select(f=frange)
    if len(data) == 0:
      raise ExperimentError('<{0}>: No data exists in range {1} - {2}'.format(self, *frange))
    return mean(data.ext)

  def adjustExtensionOffset(self, baseline, offset_range=None):
    'Adjust extension to hit baseline. If offset_range is not given, it is calculated from Pulling.extensionOffsetRange'
    offset = self.extensionOffset(offset_range) - baseline
    self.trap.ext -= offset
    self._ext_offset = offset
    return offset

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
    kwargs.setdefault('FEC', self.fits or not self.fret)
    kwargs.setdefault('title', self.filename or '')
    loc_x = min(self.trap.ext)+10
    location = list(kwargs.setdefault('annotate', (loc_x, 15)))
    if self.fret:
      self.figure.plot(self.fret, self.trap, **kwargs)
    else:
      self.figure.plot(self.trap, **kwargs)
    if self.handles:
      self.figure.plot(self.handles, hold=True)
      self.figure.annotate(unicode(self.handles), location)
      location[1] -= 1
    for fit in self.rips:
      self.figure.plot(fit, hold=True)
      text = unicode(fit)
      cutoff = 51
      if len(text) >= cutoff:
        text = text[:cutoff]+'\n'+text[cutoff:]
      self.figure.annotate(text, location)
    return self.figure

  def pickLimits(self, fig=None):
    if not fig: fig=plt.gcf()
    firstPoint,secondPoint = plt.ginput(2)

  def savefig(self, filename=None, path='.'):
    if self.figure is None or not self.figure.exists:
      raise ExperimentError('No figure available for experiment {0}'.format(str(self)))
    else: 
      filename = filename or self.filename
      self.figure.toFile(filename)

  @classmethod
  def load(cls, filename):
    return pickle.load(open(filename,'rb'))

  def save(self, filename=None):
    filename = filename or self.filename+'.exp'
    if not filename:
      raise ExperimentError('Specify a filename')
    pickle.dump( self, open(filename,'wb') )

class OpenLoop(Base):
  "Object for manipulating FretData (and optional TrapData) for open loop measurements"
  def __init__(self, pull, fret, **metadata):
    super(OpenLoop, self).__init__(pull, fret, **metadata)

  FILENAME_SYNTAX = ('min', 'sec', 'force')

  @classmethod
  def filenameMatchesType(cls, filename):
    assert isinstance(filename, str)
    finfo = fileIO.parseFilename(filename)
    if not finfo:
      logger.warning('Problem parsing filename "{}"'.format(filename))
      return False
    for attr in OpenLoop.FILENAME_SYNTAX:
      if getattr(finfo, attr) is not None:
        return True
    return False

  @classmethod
  def fromFile(cls, fretfile):
    assert isinstance(fretfile, str)
    basename, strfile, fretfile = fileIO.filesFromName(fretfile)
    if not path.exists(strfile):
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
    'filename_matching': True},
  'filtering': {
    'auto_filter': True, # apply these options automagically where needed
    'required_pulling_force': Pulling.extensionOffsetRange[1]}
  
})
