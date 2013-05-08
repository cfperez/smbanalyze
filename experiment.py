import os.path as path
from operator import itemgetter, attrgetter, methodcaller
import logging
import collections
import cPickle as pickle
import re

from numpy import where, min, max, asarray, sum, mean, all
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
  files = filter(lambda x: x.count('str')>0, fileIO.flist(*fglob))
  if not files:
    raise ExperimentError("No files found matching glob '{0}'".format(fglob))
  return fromFiles(files)

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

class List(list): #(collections.Sequence):
  def __init__(self, iterable):
    super(List, self).__init__(iterable)
    try:
      self.sort(key=attrgetter('info.mol')).sort(key=attrgetter('info.pull'))
    except AttributeError:
      pass

  def matching(self, *match):
    matched = filter(lambda x: re.search(makeMatchStrFromArgs(*match), x.filename), self)
    return List(matched)

  def get(self, name, *more):
    try:
      return map( attrgetter(name, *more), self)
    except AttributeError:
      raise ExperimentError('Missing attribute {0} in a List element'.format(action))
      
  def call(self, action, *args, **kwargs):
    try:
      return map( methodcaller(action, *args, **kwargs), self )
    except AttributeError:
      raise ExperimentError('Missing method {0} in a List element'.format(action))

  def next(self):
    try:
      return self.it.next()
    except (StopIteration, AttributeError):
      self.it = iter(self)
      return self.next()

  def plotall(self, attr, **options):
    options.setdefault('labels', self.get('filename'))
    fplot.plotall( self.get(attr), **options)

  def plot(self, **options):
    for exp in self:
      exp.plot(**options)

  def fitHandles(self, *args, **kwargs):
    return self.call('fitHandles', *args, **kwargs)

  def fitRip(self, *args, **kwargs):
    return self.call('fitRip', *args, **kwargs)

  def adjustForceOffset(self, *args, **kwargs):
    return self.call('adjustForceOffset', *args, **kwargs)

  def __getslice__(self, i, j):
    return self.__getitem__(slice(i, j))

  def __getitem__(self, key):
    out = super(List, self).__getitem__(key)
    if isinstance(out, list):
      return List(out)
    else:
      return out

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
    return plt.annotate(text, location)
    
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
    if hasattr(self,'filename'):
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

  forceOffsetRange = slice(10, 20)

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
    parameters.update(**fitOptions)
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
    mask = self.pull.maskFromLimits(x, f, (start,stop))
    fitOptions.setdefault('Lc', max(self.pull.ext)*1.05)
    fit = fitWLC(self.pull.ext, self.pull.f, mask=mask, **fitOptions)
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

  def forceOffset(self, indices=None):
    indices = indices or Pulling.forceOffsetRange
    return mean(self.pull.f[indices])

  def adjustForceOffset(self, baseline=0, offset=None):
    if offset is None:
      offset = -self.forceOffset()
    if offset != 0:
      offset -= baseline
      def geometricMean(*args):
        return 1/sum(map(lambda x: 1./x, args))
      stiffness = geometricMean(*self.metadata.get('stiffness', constants.stiffness))
      self.pull.f += offset
      self.pull.ext -= offset/stiffness
    return offset

  def recalculate(self, stiffness=None):
      if len(stiffness) != 2:
        raise ValueError('Stiffness must be 2-tuple')
      def geometricMean(*args):
        return 1/sum(map(lambda x: 1./x, args))
      current_k = self.metadata.get('stiffness', constants.stiffness)
      new_k = stiffness or current_k
      mean_k = geometricMean(*new_k)
      beadRadii = self.metadata.get('bead_radii', constants.sumOfBeadRadii)
      self.pull.f *= min(new_k)/min(current_k)
      self.pull.ext = self.pull.sep - beadRadii - self.pull.f/mean_k

  def extensionOffset(self, frange=None):
    'Returns average extension of FEC between given forces'
    frange = frange or Pulling.frange

    return mean(self.pull.ext[ext])

  def adjustExtensionOffset(self, baseline, offset=None):
    pass
    
  def plot(self, **kwargs):
    kwargs.setdefault('FEC', self.fits or not self.fret)
    kwargs.setdefault('title', self.filename or '')
    loc_x = min(self.pull.ext)+10
    kwargs.setdefault('location', (loc_x, 15))
    self.figure.plot(self.fret, self.pull, **kwargs)
    if self.handles:
      self.figure.plot(self.handles, hold=True)
    for fit in self.rips:
      self.figure.plot(fit, hold=True)
      text = str(fit)
      self.figure.annotate(text, kwargs['location'])

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
