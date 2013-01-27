import os
import operator
import glob
import logging

import numpy as np

import FileIO 
import FRET
import Image
import Constants
from Types import *

class ExperimentError(Exception):
  pass

def fromFiles(*filelist):
  return [Pulling.fromFile(fname) for fname in filelist]

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, pull speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  def __init__(self, data, **metadata):
    self._data=data
    self.metadata=metadata
    self.hasfret = hasFretData(self._data)

  @property
  def _fields(self):
    return self._data._fields

  def __getattr__(self,attr):
    if attr in self._fields:
      return getattr(self._data,attr)
    raise AttributeError("'%s' has no attribute %s" % 
      (self.__class__.__name__,attr))

  def plot(self):
    raise NotImplementedError

  def plot(self):
    raise NotImplementedError

  @property
  def molname(self):
    return FRET.molname(self)

class Pulling(Base):
  "stretching curves of single molecule .str camera .cam image .img"
  # pull = Pulling(...)
  # pull.molname = 's1m4'
  # pull.construct = 'SJF'
  # pull.conditions = '1B 1.0 nM'
  #
  # pull.plot(**kwargs) = sets up good defaults for title, names, etc.
  # pull.ext pull.f pull.sep pull.fret pull.donor

  def __init__(self,pull,fret=None, **metadata):
    if hasPullFretData(pull) or fret is None:
      super(Pulling,self).__init__(pull, **metadata)
    elif fret is not None:
      super(Pulling,self).__init__(PullFretData._make(pull+fret), **metadata)


  @classmethod
  def fromFile(cls,strfile,fretfile=None):
    basename,ext=os.path.splitext(strfile)
    if not ext:
      strfile = FileIO.add_pull_ext(basename)

    # check if base + .fret exists if not specified already
    # and use it, or else load/don't load fretfile
    fretfileFromBase = FileIO.add_fret_ext(basename)
    if not fretfile and os.path.exists(fretfileFromBase):
      fret = FileIO.load(fretfileFromBase)
    else:
      fret = fretfile and FileIO.load(fretfile)
    meta,data = FileIO.load(strfile,commentparser=FileIO.commentsToSettings)
    newPull = cls(data,fret,**meta)
    newPull.file = basename
    try:
      newPull.info = FileIO.parseFilename(basename)
    except:
      logging.warning('Problem parsing filename %s' % basename)

    return newPull

  def plot(self, **kwargs):
    kwargs.setdefault('FEC',not self.hasfret)
    FRET.plot(self._data, **kwargs)

  def __len__(self):
    return self.size

  def __iter__(self):
    return iter(self._array)

  @property
  def molname(self):
    return 's{}m{}'.format(self.slide,self.mol)

class OpenLoop(Base):
  "camera .cam image. img"
  pass

class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass
