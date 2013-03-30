import os
from operator import itemgetter
import logging
import collections

from numpy import where, min, max
import numpy as np
import matplotlib.pyplot as plt

import fileIO 
from curvefit import Fit
from useful import isInt, groupat
import fplot
import constants
from types import *

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
    for pull,fret in groupat(hasPullData, datalist, size=2):
      output += [Pulling(pull,fret)]
    return output if len(output)>1 else output[-1]

def fromMatch(*fglob):
  files = filter(lambda x: x.count('str')>0,fileIO.flist(*fglob))
  if not files:
    raise ExperimentError("No files found matching glob '{0}'".format(fglob))
  return fromFiles(files)

def fromFiles(*filelist):
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
      base, ext = os.path.splitext(filename)
      if ext[1:] not in ('emf', 'eps', 'pdf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz'):
        filename += constants.DEFAULT_FIGURE_EXT
    else:
      filename = 'Figure {0}{1}'.format(self.figure.number, constants.DEFAULT_FIGURE_EXT)

    self.figure.savefig(filename, bbox_inches='tight', pad_inches=0.1)

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, pull speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  def __init__(self, data, **metadata):
    self._data=data
    self.metadata=metadata
    self.hasfret = hasFretData(self._data)

  @property
  def fields(self):
    return self._data._fields

  def __getattr__(self,attr):
    if attr in self.fields:
      return getattr(self._data,attr)
    raise AttributeError("'%s' has no attribute %s" % 
      (self,attr))

  def __repr__(self):
    if hasattr(self,'filename'):
      return "<Experiment.%s from '%s'>" % (self.__class__.__name__, self.filename)
    else:
      return super(Base,self).__repr__()

  @property
  def pullData(self):
    return PullData(*self._data[0:3])

  @property
  def fec(self):
    try:
      return np.asarray((self._data.ext,self._data.f)).T
    except AttributeError:
      raise ExperimentError('Experiment %s does not have FEC data' % self)

  def plot(self):
    raise NotImplementedError
    
class Pulling(Base):
  "stretching curves of single molecule .str camera .cam image .img"

  def __init__(self,pull,fret=None, **metadata):
    if hasPullFretData(pull) or fret is None:
      super(Pulling,self).__init__(pull, **metadata)
    elif fret is not None:
      super(Pulling,self).__init__(PullFretData._make(pull+fret), **metadata)

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
    meta,data = fileIO.load(strfile,comments=fileIO.toSettings)

    # check if base + .fret exists if not specified already
    # and use it, or else load/don't load fretfile
    fretfileFromBase = fileIO.add_fret_ext(basename)
    if not fretfile and os.path.exists(fretfileFromBase):
      cam_metadata, fret = fileIO.load(fretfileFromBase, comments=fileIO.toSettings)
      meta.update(cam_metadata)
    elif fretfile and os.path.exists(fretfile):
      cam_metadata, fret = fretfile and fileIO.load(fretfile, comments=fileIO.toSettings)
      meta.update(cam_metadata)
    else: fret=None

    newPull = cls(data,fret,**meta)
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
    fitOptions.setdefault('fitfunc', 'MMS')

    ext_fit, f_fit = self._constrainFitDataFromLimits(x, f, (start,stop))

    if len(ext_fit)==0 or len(f_fit)==0:
      raise ExperimentError('Provided constraints (%s, %s) are outside of data' %
            (x, f))

    fitOptions.setdefault('Lc', max(ext_fit)*1.05)
    try:
      if fitOptions['fitfunc'].startswith('MMS'):
        init = (fitOptions.get('fixed', ()))
        fitOptions['fixed']=  init+('K',)
    except AttributeError: pass

    fit = Fit(ext_fit, f_fit, **fitOptions)
    self.lastFit = fit

    if self.figure.exists:
      self.figure.plot(fit, hold=True)

    return fit

  def _constrainFitDataFromLimits(self, x, f, limits=(0,-1)):
    start, stop = limits
    ext_fit,f_fit = self.ext[start:stop], self.f[start:stop]

    if f is None: f=[np.max(f_fit)]
    try:
      min_f, max_f = f
    except TypeError:
      min_f, max_f = min(f_fit), f
    except ValueError:
      min_f, max_f = min(f_fit), f[0]
    if max_f>max(f_fit): max_f=max(f_fit)

    if x is None: x=[min(ext_fit)]
    try:
      min_ext, max_ext = x
    except TypeError:
      min_ext, max_ext = x, ext_fit[min(where(f_fit>=max_f)[0])-1]
    except ValueError:
      min_ext, max_ext = x[0], ext_fit[min(where(f_fit>=max_f)[0])-1]

    between = lambda s,a,b: (s>=a) & (s<=b)
    mask = between(ext_fit, min_ext, max_ext)

    return ext_fit[mask], f_fit[mask]

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
        return 1/np.sum(map(lambda x: 1./x, args))
      beadRadii = self.metadata.get('bead_radii', constants.sumOfBeadRadii)
      stiffness = geometricMean(*self.metadata.get('stiffness', constants.stiffness))
      self.f -= offset
      self.ext = self.sep - beadRadii - self.f/stiffness

  def plot(self, **kwargs):
    kwargs.setdefault('FEC', self.fits or not self.hasfret)
    title=self.filename or ''
    self.figure.plot(self._data, title=title, **kwargs)
    for fit in self.fits:
      self.figure.plot(fit, hold=True)

  def pickLimits(fig=None):
    if not fig: fig=plt.gcf()
    firstPoint,secondPoint = ginput(2)

  def savefig(self, filename=None, path='.'):
    if not self.figure or not self.figure.exists:
      raise ExperimentError('No figure available for object {0}'.format(self))
    else: 
      filename = filename or self.filename
      self.figure.toFile(filename)

class OpenLoop(Base):
  "camera .cam image. img"
  pass

class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass
