import os.path as path
from operator import itemgetter
import logging
import collections

from numpy import where, min, max, asarray, sum
import matplotlib.pyplot as plt

import fileIO 
from curvefit import Fit, fitWLC
from useful import isInt, groupat
import fplot
import constants
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
  files = filter(lambda x: x.count('str')>0,fileIO.flist(*fglob))
  if not files:
    raise ExperimentError("No files found matching glob '{0}'".format(fglob))
  return fromFiles(files)

def fromFiles(*filelist):
  'Load experiments from a list of files or an argument list'
  filelist = collapseArgList(filelist)
  return [Pulling.fromFile(fname) for fname in filelist]

def collapseArgList(arglist):
  'Allows a function to take either *args or a list as its argument'
  if isinstance(arglist[0], list) or isinstance(arglist[0], tuple):
    if len(arglist)>1: raise ValueError('Argument must be a list or *args')
    return arglist[0]
  else:
    return arglist

class ExperimentList(collections.Sequence):
  def __init__(self):
    pass

  def __len__(self):
    pass

  def __iter__(self):
    pass

  def __contains__(self):
    pass

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
    self.figure.clf()
    self.figure.show()

  def toFile(self, filename=None):
    if filename:
      base, ext = path.splitext(filename)
      if ext[1:] not in ('emf', 'eps', 'pdf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz'):
        filename += constants.DEFAULT_FIGURE_EXT
    else:
      filename = 'Figure {0}{1}'.format(self.figure.number, constants.DEFAULT_FIGURE_EXT)

    self.figure.savefig(filename, bbox_inches='tight', pad_inches=0.1)

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, pull speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  def __init__(self, pull, fret, **metadata):
    self.fec = pull
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
    
class Pulling(Base):
  "stretching curves of single molecule .str camera .cam image .img"

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
    parameters.update(Lc1=guess)
    parameters.setdefault('K', 1100)
    fitOptions.update( {'fixed': ('K','K1','Lc','Lp','F0','Lp1'), 
          'fitfunc': 'MMS_rip'}, **parameters)
    rip = self.fitForceExtension(x, f, **fitOptions)
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
    mask = self.fec.maskFromLimits(x, f, (start,stop))
    fitOptions.setdefault('Lc', max(self.fec.ext)*1.05)
    fit = fitWLC(self.fec.ext, self.fec.f, mask=mask, **fitOptions)
    self.lastFit = fit

    if self.figure.exists:
      self.figure.plot(fit, hold=True)

    return fit

  def fitFEC(self, x=None, f=None, tolerance=0.5, **fitOptions):
    offset = 0
    loops = 0
    while True:
      loops += 1
      self.adjustForceOffset(offset)
      fit=self.fitForceExtension(x,f,**fitOptions)
      offset = abs(fit['F0'])
      if offset < tolerance or loops>10: break
    return fit
    
  @property
  def fits(self):
    if self.handles:
      return [self.handles]+self.rips
    else:
      return [self.lastFit] if self.lastFit else []

  def adjustForceOffset(self, offset):
    if offset!=0:
      def geometricMean(*args):
        return 1/sum(map(lambda x: 1./x, args))
      beadRadii = self.metadata.get('bead_radii', constants.sumOfBeadRadii)
      stiffness = geometricMean(*self.metadata.get('stiffness', constants.stiffness))
      self.f -= offset
      self.ext = self.sep - beadRadii - self.f/stiffness

  def plot(self, **kwargs):
    kwargs.setdefault('FEC', self.fits or not self.fret)
    kwargs.setdefault('title', self.filename or '')
    self.figure.plot(self.fret, self.fec, **kwargs)
    for fit in self.fits:
      self.figure.plot(fit, hold=True)

  def pickLimits(fig=None):
    if not fig: fig=plt.gcf()
    firstPoint,secondPoint = ginput(2)

  def savefig(self, filename=None, path='.'):
    if self.figure is None or not self.figure.exists:
      raise ExperimentError('No figure available for experiment {0}'.format(str(self)))
    else: 
      filename = filename or self.filename
      self.figure.toFile(filename)

class OpenLoop(Base):
  "camera .cam image. img"
  pass

class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass
