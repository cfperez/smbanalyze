import pdb
import os
import operator
import glob
import logging
import collections

from numpy import where, min, max
import numpy as np
import matplotlib.pyplot as plt

import FileIO 
from curvefit import Fit
from useful import isInt, groupat
import FRET
import Image
import Constants
from Types import *

logger = logging.getLogger(__name__)
logger.setLevel(Constants.logLevel)
logger.addHandler(Constants.logHandler)

class ExperimentError(Exception):
  pass

def fromGlob(*globs, **kwargs):
  exptype = kwargs.get('type','pull')
  if exptype=='refold':
    filenames = FileIO.flist(globs)

# Good abstraction for above
#Experiment.fromData(pull_data, [fret_data,] type='pull')
#Experiment.fromData(pull_fret_data_1, pull_fret_data_2, type='loopdelay')
# RETURNS Experiment subclass that can:
# 1) Keep track and make sense of metadata (loop times, stiffnesses)
# 2) Provide structure/context for analysis functions that need that data
# 3) Convenience functions for plotting and saving
# 4) Ability to stack with other Experiments for global analysis
# 5) Enable potential to develop summary outputs that save to database/file (external functions)

def fromData(*datalist, **kwargs):
  "List of experiments from PullData and FretData type"
  exptype = kwargs.get('type','pull')
  if exptype=='pull':
    output = []
    for pull,fret in groupat(hasPullData, datalist, size=2):
      output += [Pulling(pull,fret)]
    return output if len(output)>1 else output[-1]

def fromMatch(*fglob):
  return fromFiles(filter(lambda x: x.count('str')>0,FileIO.flist(*fglob)))

def fromFiles(filelist):
  return [Pulling.fromFile(fname) for fname in filelist]

class ExperimentList(collections.Sequence):
  def __init__(self):
    pass

  def __len__(self):
    pass

  def __iter__(self):
    pass

  def __contains__(self):
    pass

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
    if hasattr(self,'file'):
      return "<Experiment.%s from '%s'>" % (self.__class__.__name__, self.file)
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

    self.figure = None
    self.fit = None

  @classmethod
  def fromFile(cls,strfile,fretfile=None):
    basename,ext=FileIO.splitext(strfile)
    if not ext:
      strfile = FileIO.add_pull_ext(basename)
    meta,data = FileIO.load(strfile,comments=FileIO.toSettings)

    # check if base + .fret exists if not specified already
    # and use it, or else load/don't load fretfile
    fretfileFromBase = FileIO.add_fret_ext(basename)
    if not fretfile and os.path.exists(fretfileFromBase):
      cam_metadata, fret = FileIO.load(fretfileFromBase, comments=FileIO.toSettings)
      meta.update(cam_metadata)
    elif fretfile and os.path.exists(fretfile):
      cam_metadata, fret = fretfile and FileIO.load(fretfile, comments=FileIO.toSettings)
      meta.update(cam_metadata)
    else: fret=None

    newPull = cls(data,fret,**meta)
    newPull.file = basename
    try:
      newPull.info = FileIO.parseFilename(basename)
    except:
      logger.warning('Problem parsing filename %s' % basename)

    return newPull

  def plot(self, **kwargs):
    kwargs.setdefault('FEC', not self.hasfret)
    FRET.plot(self._data, **kwargs)
    self.figure = plt.gcf()
    if self.fit:
      self.fit.plot(hold=True)

  def fitRip(self, x=None, f=None, start=0, stop=-1, **fitOptions):
    pass
    
  def fitForceExtension(self, x=None, f=None, start=0, stop=-1, **fitOptions):
    fitOptions.setdefault('fitfunc', 'MS')

    ext_fit, f_fit = self._constrainFitDataFromLimits(x, f, (start,stop))

    fitOptions.setdefault('Lc', max(ext_fit)*1.05)
    try:
      if fitOptions['fitfunc'].startswith('MMS'):
        fitOptions.setdefault('fixed', 'K')
    except AttributeError: pass

    self.fit = Fit(ext_fit, f_fit, **fitOptions)

    if self.figure and plt.fignum_exists(self.figure.number):
      plt.figure(self.figure.number)
      self.fit.plot(hold=True)

    return self.fit

  def _constrainFitDataFromLimits(self, x, f, limits):
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

    if x is None or x<min(ext_fit): x=[min(ext_fit)]
    try:
      min_ext, max_ext = x
    except TypeError:
      min_ext, max_ext = x, ext_fit[min(where(f_fit>=max_f)[0])-1]
    except ValueError:
      min_ext, max_ext = x[0], ext_fit[min(where(f_fit>=max_f)[0])-1]

    between = lambda s,a,b: (s>=a) & (s<=b)
    mask = between(ext_fit, min_ext, max_ext)

    return ext_fit[mask], f_fit[mask]

  def adjustForceOffset(self,offset):
    def geometricMean(*args):
      return 1/np.sum(map(lambda x: 1./x, args))
    beadRadii = self.metadata.get('bead_radii', Constants.sumOfBeadRadii)
    stiffness = geometricMean(*self.metadata.get('stiffness', Constants.stiffness))
    self.f -= offset
    self.ext = self.sep - beadRadii - self.f/stiffness

  def pickLimits(fig=None):
    if not fig: fig=plt.gcf()
    firstPoint,secondPoint = ginput()

class OpenLoop(Base):
  "camera .cam image. img"
  pass

class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass
