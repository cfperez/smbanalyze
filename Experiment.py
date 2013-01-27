import os
import operator
import glob

import numpy as np

import FileIO
import FRET
import Image
import Constants
from Types import *

class ExperimentError(Exception):
  pass

def loadPull(fileglob, **kwargs):
  "Create a new Experiment subclass based on filetype and data"
  verbose = kwargs.get('verbose')
  recalc = kwargs.get('recalc')
  roi = kwargs.get('roi_file')
  #autoglob = kwargs.get('autoglob', True)
  exact_match = kwargs.get('exact_match', False)

  fret = []
  pulldata = []
  finfo = []
  bg_files = None

  if os.path.isfile(fileglob):
    files = [fileglob]
  else:
    files = (exact_match and glob.glob(fileglob)) or glob.glob('*%s[._]*str'%fileglob)
  if not files:
    raise IOError("No files found using glob '%s'" % fileglob)

  for filename in files:

    basename,ext = os.path.splitext(filename)
    fretfile = FileIO.add_fret_ext(basename)
    imgfile = FileIO.add_img_ext(basename)

    if recalc and not fret:
      print "-- Recalculating FRET data --"
      FRET.calcDirectory(fileglob, **kwargs)

    if verbose:
      print "Experiment.load: Processing %s..." % basename

    if os.path.isfile(fretfile):
      if verbose:
        print "\tLoading fret from " + fretfile
      fret += [FileIO.loadfret(fretfile)]
    elif os.path.isfile(imgfile):
      if verbose:
        print "\tCalculating fret from " + imgfile

      if not bg_files:
        matched = FRET.matchImgFilestoBackground()
        bg_files = dict([ (img,bg) for img,bg in matched if fileglob in img ])
        if roi:
          roi = Image.ROI.fromFile(roi)
          Image.setDefaultROI(*roi)

      image = Image.Stack(imgfile)
      image -= Image.fromBackground(bg_files[imgfile])
      fret += [(image.time,image.donor,image.acceptor,
                FRET.calcToFile(image,fretfile))]
    else:
      if verbose:
        print "\tNo .fret or .img file found"
      fret += [None]

    pullfile = FileIO.add_pull_ext(basename)
    if verbose:
      print "\tLoading pulling data from %s" % pullfile
    pulldata += [FileIO.loadstr(pullfile)]
    
    finfo += [FileIO.parseFilename(filename)]

  return Pulling(finfo,pulldata,fret)

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, pull speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  def __init__(self, data, fields=None):
    self._data=data
    self._array=np.asarray(self._data).T

  @property
  def _fields(self):
    return self._data._fields

  def __getattr__(self,attr):
    if attr in self._fields:
      return getattr(self._data,attr)
    raise AttributeError("'%s' has no attribute %s" % 
      (self.__class__.__name__,attr))

  def _setattr(self,named):
    for attr,val in named._asdict().iteritems():
      self.fieldnames += (attr,)
      setattr(self, attr, val)

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

  def __init__(self,pull,fret=None):
    if hasPullFretData(pull) or fret is None:
      super(Pulling,self).__init__(pull)
    elif fret is not None:
      super(Pulling,self).__init__(PullFretData._make(pull+fret))

    self.hasfret = hasFretData(self._data)

  def plot(self, **kwargs):
    FRET.plot(self._data, **kwargs)

  def __getitem__(self,key):
    return self._array[key]

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
