import os.path as path
from operator import itemgetter, attrgetter, methodcaller, add
import logging
import cPickle as pickle
import re

from numpy import min, max, asarray, sum, mean, all, any
import matplotlib.pyplot as plt

import fileIO 
from curvefit import fitWLC
from useful import groupat, makeMatchStrFromArgs
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

def fromMatch(*fglob):
  'Load experiments from files using concatenation of argument list as a glob'
  fglob = list(fglob)
  files = fileIO.flist(fglob[0], fileIO.PULL_FILE)
  fglob[0] = path.basename(fglob[0])
  if not files:
    raise ExperimentError("No files found matching glob '{0}'".format(fglob))
  matched = filter(lambda x: re.search(makeMatchStrFromArgs(*fglob), fileIO.splitext(x)[0]), files)
  return fromFiles(matched)

def fromFiles(*filelist):
  'Load experiments from a list of files or an argument list'
  filelist = collapseArgList(filelist)
  return List(map(fromFile, filelist))

def fromFile(filename):
  if OpenLoop.hasFiletype(filename):
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
  def __init__(self, iterable):
    super(List, self).__init__(iterable)
    try:
      self.sort(key=attrgetter('info.mol'))
      self.sort(key=attrgetter('info.pull'))
    except AttributeError:
      pass

  def filter(self, condition):
    "Returns a List with experiments matching condition"
    return List( filter(condition, self) )

  def matching(self, *match):
    "Return List with experiment names matching *match globs"
    return self.filter(lambda x: re.search(makeMatchStrFromArgs(*match), x.filename))

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
    self.figure = Figure.fromCurrent()
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

  def fitRip(self, *args, **kwargs):
    return self.call('fitRip', *args, **kwargs)

  def adjustForceOffset(self, *args, **kwargs):
    return self.call('adjustForceOffset', *args, **kwargs)

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
    foffset = self.adjustForceOffset(to_f, x_range)
    mean_x = to_x or mean(self.call('extensionOffset', f_range))
    xoffset = self.adjustExtensionOffset(mean_x, f_range)
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

class Figure(object):
  def __init__(self, fignumber=None):
    self.figure = plt.figure(fignumber) if fignumber else None

  @classmethod
  def fromCurrent(cls):
    return cls(plt.gcf().number)

  def new(self):
    self.figure = plt.figure()

  @property
  def exists(self):
    return self.figure is not None and plt.fignum_exists(self.figure.number)

  def makeCurrent(self):
    if not self.exists:
      raise RuntimeError('Figure object does not exist')
    plt.figure(self.figure.number)
    return self

  def plot(self, *args, **kwargs):
    if not self.exists:
      self.new()
    else:
      self.makeCurrent()
    try:
      # Treat the first argument as an object that can plot itself...
      return args[0].plot(*args[1:], **kwargs)
    except AttributeError:
      # ...unless it can't
      return fplot.plot(*args, **kwargs)

  def clear(self):
    if self.exists:
      self.figure.clf()
      self.figure.show()

  def annotate(self, text, location):
    "Annotate figure with text at location (x,y)"
    x,y = location
    return plt.text(x, y, text)
  IMAGE_OUTPUT_FORMATS = ('emf', 'eps', 'pdf', 'png', 'ps',
      'raw', 'rgba', 'svg', 'svgz') 

  DEFAULT_SIZE = (9, 7.5)
  def toFile(self, filename=None):
    if filename:
      ext = path.splitext(filename)[1]
      if ext[1:] not in Figure.IMAGE_OUTPUT_FORMATS:
        filename += constants.DEFAULT_FIGURE_EXT
    else:
      filename = 'Figure {0}{1}'.format(self.figure.number, constants.DEFAULT_FIGURE_EXT)
    self.figure.set_size_inches(*Figure.DEFAULT_SIZE)
    self.figure.savefig(filename, bbox_inches='tight', pad_inches=0.1)

class TrapSettings(object):
  pass

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, trap speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
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
    self.figure = Figure()
    self.filename = ''

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
    self.rips = []
    self.lastFit = None
    self._appendNewRip = True
    self._ext_offset = 0

  @classmethod
  def aggregate(cls, *pulling):
    pass

  @classmethod
  def fromFile(cls, strfile, fretfile=None):
    'Load stretching data and corresponding fret data from files'
    basename, strfile, fretfileFromBase = fileIO.filesFromName(strfile)
    fretfile = fretfile or fretfileFromBase
    if not path.exists(fretfile):
      fretfile = None
    return super(Pulling, cls).fromFile(strfile, fretfile)

  def loadimg(self, path='.', **kwargs):
    filename = self.filename
    try:
      return image.fromFile(filename, **kwargs)
    except IOError:
      raise ExperimentError('IOError loading file: check image file location!')

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

  @property
  def ripSizes(self):
    rips = asarray(map(itemgetter('Lc1'), self.rips))
    if len(rips) == 1:
      return rips
    else:
      return rips[1:]-rips[:-1]

  def fitForceExtension(self, x=None, f=None, start=0, stop=-1, **fitOptions):
    "Fit WLC to Pulling curve and plot"
    pull = self.trap.select(x, f, (start,stop))
    fit = fitWLC(pull.ext, pull.f, **fitOptions)
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
      raise ExperimentError('No data exists in range {0} - {1}'.format(*x_range))
    return mean(data.f)

  @property
  def meanStiffness(self, stiffness=None):
    inverseAverage = lambda args: 1/sum(map(lambda x: 1./x, args))
    return inverseAverage(self.trap.metadata.get('stiffness', constants.stiffness))

  def adjustForceOffset(self, baseline=0.0, offset=None, offset_range=None):
    offset = offset or -self.forceOffset(offset_range)
    offset += baseline
    self.trap.f += offset
    self.trap.ext -= offset/self.meanStiffness
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

def hasAnyAttr(obj, *attr):
  return any( map(lambda a: hasattr(obj, a), attr) )

class OpenLoop(Base):
  "Object for manipulating FretData (and optional TrapData) for open loop measurements"
  def __init__(self, pull, fret, **metadata):
    super(OpenLoop, self).__init__(pull, fret, **metadata)

  OPENLOOP_FNAME_SYNTAX = ('min', 'sec', 'force')

  @classmethod
  def hasFiletype(cls, filename):
    assert isinstance(filename, str)
    finfo = fileIO.parseFilename(filename)
    if not finfo:
      logger.warning('Problem parsing filename "{}"'.format(filename))
      return False
    for attr in OpenLoop.OPENLOOP_FNAME_SYNTAX:
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
