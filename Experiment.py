import os
import operator
import glob
import logging
import collections

import numpy as np
import matplotlib.pyplot as plt

import FileIO 
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

# Harder for loop delay experiments
# >>> Experiment.fromFiles('s1m1.img','s1m1.str','s1m1_2_refold_1.0s.img',
# 's1m1_2_up.img','s1m1_2_up.str')

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
  return fromFiles(FileIO.flist(*fglob))

def fromFiles(*filelist):
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

  @classmethod
  def fromFile(cls,strfile,fretfile=None):
    basename,ext=FileIO.splitext(strfile)
    if not ext:
      strfile = FileIO.add_pull_ext(basename)

    # check if base + .fret exists if not specified already
    # and use it, or else load/don't load fretfile
    fretfileFromBase = FileIO.add_fret_ext(basename)
    camfileFromBase = FileIO.add_cam_ext(basename)
    if not fretfile and os.path.exists(fretfileFromBase):
      fret = FileIO.load(fretfileFromBase)
    else:
      fret = fretfile and FileIO.load(fretfile)
    meta,data = FileIO.load(strfile,comments=FileIO.toSettings)
    meta.update(FileIO.loadcam(camfileFromBase))

    newPull = cls(data,fret,**meta)
    newPull.file = basename
    try:
      newPull.info = FileIO.parseFilename(basename)
    except:
      logger.warning('Problem parsing filename %s' % basename)

    return newPull

  def plot(self, **kwargs):
    kwargs.setdefault('FEC',not self.hasfret)
    FRET.plot(self._data, **kwargs)

  def pickLimits(fig=None):
    if not fig: fig=plt.gcf()
    firstPoint,secondPoint = ginput()

class OpenLoop(Base):
  "camera .cam image. img"
  pass

class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass
