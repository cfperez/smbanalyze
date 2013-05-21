import os.path as path
from operator import itemgetter, attrgetter, methodcaller
import logging
import collections
import cPickle as pickle
import re

from numpy import where, min, max, asarray, sum, mean, all, linspace
import matplotlib.pyplot as plt

import fileIO 
from curvefit import Fit, fitWLC
from useful import isInt, groupat, makeMatchStrFromArgs
import fplot
import constants
import image
from datatypes import *

logger = logging.getLogger(__name__)
logger.setLevel(constants.logLevel)
logger.addHandler(constants.logHandler)

class ExperimentError(Exception):
  pass

def fromData(*datalist, **kwargs):
  "List of experiments from PullData and FretData type"
  exptype = kwargs.get('type','pull')
  if exptype=='pull':
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
  basenames = map(lambda x: fileIO.splitext(x)[0], files)
  matched = filter(lambda x: re.search(makeMatchStrFromArgs(*fglob), fileIO.splitext(x)[0]), files)
  return fromFiles(matched)

def fromFiles(*filelist):
  'Load experiments from a list of files or an argument list'
  filelist = collapseArgList(filelist)
  return List([Pulling.fromFile(fname) for fname in filelist])

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
    "Get all attributes of experiments in List by name: get('pull','f') => pull.f"
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
      
  def call(self, action, *args, **kwargs):
    "Call function <action> with *args and **kwargs on all experiments"
    try:
      return map( methodcaller(action, *args, **kwargs), self )
    except AttributeError:
      raise ExperimentError('Missing method {0} in a List element'.format(action))

  METADATA_CHECK = {'pull': ('step_size', 'sampling_time'),
                    'fret': ('exposurems', 'frames', 'gain', 'binning') }

  def aggregate(self, datatype='pull', check_metadata=None):
    "Aggregate data of given kind from experiments by appending."
    if check_metadata is not None:
      assert isinstance(check_metadata, tuple) or isinstance(check_metadata, list)

    filtered = self.has(datatype)
    logger.info('Using experiments {}'.format(filtered))
    aggregated = sum(filtered.get(datatype)) or None

    meta = filtered.get('metadata')
    check_for_fields = check_metadata or List.METADATA_CHECK[datatype]
    for field in check_for_fields:
      try:
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
    fec, fdata = None, None
    fec = self.aggregate('pull')
    fdata = self.aggregate('fret')
    return Pulling(fec, fdata)

  def plotall(self, attr=None, **options):
    "Plot all experiments overlayed onto same panel. Can plot only pull or fret."
    assert attr in set([None,'pull','fret'])
    options.setdefault('labels', self.get('filename'))
    if hasattr(self,'figure') and self.figure.exists:
     self.figure.makeCurrent()
     self.figure.clear()
    if attr is None and self.has('fret'):
      options.setdefault('FEC', False)
      options.setdefault('legend', None)
      fplot.plotall(self.get('fret'), self.get('pull'), hold=True, **options)
    else:
      attr = attr or 'pull'
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
    
  def toFile(self, filename=None):
    if filename:
      base, ext = path.splitext(filename)
      if ext[1:] not in ('emf', 'eps', 'pdf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz'):
        filename += constants.DEFAULT_FIGURE_EXT
    else:
      filename = 'Figure {0}{1}'.format(self.figure.number, constants.DEFAULT_FIGURE_EXT)

    self.figure.set_size_inches(9,7.5)
    self.figure.savefig(filename, bbox_inches='tight', pad_inches=0.1)

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, pull speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  def __init__(self, pull, fret, **metadata):
    self.pull = pull
    self.fret = fret
    self.metadata = pull.metadata
    if fret:
      self.metadata.update(fret.metadata)
    self.metadata.update(metadata)

  def __repr__(self):
    if getattr(self,'filename', None):
      return "<Experiment.%s from '%s'>" % (self.__class__.__name__, self.filename)
    else:
      return super(Base,self).__repr__()

  def __str__(self):
    return self.__repr__()

  def plot(self):
    raise NotImplementedError
    
def add_to_all(arg, x):
  try:
    return arg+x
  except TypeError:
    return map(lambda y: y+x, arg)

class Pulling(Base):
  "stretching curves of single molecule .str camera .cam image .img"

  forceOffsetRange = (740,770)
  extensionOffsetRange = (13,16)

  def __init__(self, pull, fret=None, **metadata):
    if not hasTrapData(pull):
      raise ExperimentError("Argument 'pull' must contain trapping data")
    if fret and not hasFretData(fret):
      raise ExperimentError("Argument 'fret' must contain FRET data")

    super(Pulling, self).__init__(pull, fret, **metadata)

    self.figure = Figure()
    self.handles = None
    self.rips = []
    self.lastFit = None
    self._appendNewRip = True
    self.filename = ''
    self._ext_offset = 0

  @classmethod
  def fromFile(cls,strfile,fretfile=None):
    'Load stretching data and corresponding fret data from files'
    basename,ext=fileIO.splitext(strfile)
    if not ext:
      strfile = fileIO.add_pull_ext(basename)
    pull = TrapData.fromFile(strfile)
    metadata = {}

    # check if base + .fret exists if not specified already
    # and use it, or else load/don't load fretfile
    fretfileFromBase = fileIO.add_fret_ext(basename)
    if not fretfile and path.exists(fretfileFromBase):
      fretfile = fretfileFromBase
    elif fretfile and not path.exists(fretfile):
      raise ExperimentError("Fret file {0} not found".format(fretfile))
    fret = FretData.fromFile(fretfile) if fretfile else None

    newPull = cls(pull, fret, **metadata)
    newPull.filename = basename

    try:
      newPull.info = fileIO.parseFilename(basename)
    except:
      logger.warning('Problem parsing filename %s' % basename)

    return newPull

  def loadimg(self, path='.', **kwargs):
    filename = self.filename
    return image.Stack(fileIO.add_img_ext(filename))
    return image.fromFile(filename, **kwargs)

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
    return map(itemgetter('Lc1'), self.rips)

  def fitForceExtension(self, x=None, f=None, start=0, stop=-1, **fitOptions):
    "Fit WLC to Pulling curve and plot"
    pull = self.pull.select(x, f, (start,stop))
    fitOptions.setdefault('Lc', max(self.pull.ext)*1.05)
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

  def forceOffset(self, xrange=None):
    xrange = xrange or Pulling.forceOffsetRange
    data = self.pull.select(x=xrange)
    if len(data) == 0:
      raise ExperimentError('No data exists in range {0} - {1}'.format(*xrange))
    return mean(data.f)

  def adjustForceOffset(self, baseline=0.0, offset=None, offset_range=None):
    offset = offset or -self.forceOffset(offset_range)
    offset += baseline
    def geometricMean(*args):
      return 1/sum(map(lambda x: 1./x, args))
    stiffness = geometricMean(*self.metadata.get('stiffness', constants.stiffness))
    self.pull.f += offset
    self.pull.ext -= offset/stiffness
    return offset

  def extensionOffset(self, frange=None):
    'Returns average extension of FEC between given forces'
    frange = frange or Pulling.extensionOffsetRange
    data = self.pull.select(f=frange)
    if len(data) == 0:
      raise ExperimentError('<{0}>: No data exists in range {1} - {2}'.format(self, *frange))
    return mean(data.ext)

  def adjustExtensionOffset(self, baseline, offset_range=None):
    'Adjust extension to hit baseline. If offset_range is not given, it is calculated from Pulling.extensionOffsetRange'
    offset = self.extensionOffset(offset_range) - baseline
    self.pull.ext -= offset
    self._ext_offset = offset
    return offset

  def recalculate(self, stiffness=None):
      if stiffness and len(stiffness) != 2:
        raise ValueError('Stiffness must be 2-tuple')
      current_k = self.metadata.get('stiffness', constants.stiffness)
      ratio_current_k = min(current_k)/max(current_k)
      new_k = stiffness or current_k
      self.metadata['stiffness'] = new_k
      beadRadii = self.metadata.get('bead_radii', constants.sumOfBeadRadii)

      displacement = self.pull.f/min(current_k)
      ratio = 1+min(new_k)/max(new_k)

      self.pull.f = displacement*min(new_k)
      self.pull.ext = self.pull.sep - beadRadii - displacement*ratio - self._ext_offset
      return self

  def plot(self, **kwargs):
    kwargs.setdefault('FEC', self.fits or not self.fret)
    kwargs.setdefault('title', self.filename or '')
    loc_x = min(self.pull.ext)+10
    location = list(kwargs.setdefault('annotate', (loc_x, 15)))
    if self.fret:
      self.figure.plot(self.fret, self.pull, **kwargs)
    else:
      self.figure.plot(self.pull, **kwargs)
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

  def pickLimits(fig=None):
    if not fig: fig=plt.gcf()
    firstPoint,secondPoint = ginput(2)

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
  "camera .cam image. img"
  pass

class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass
