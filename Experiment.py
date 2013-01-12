import os
import operator
import glob

import numpy as np

from useful import dotdict
import FileIO
import FRET
import Image
import Constants
from Types import FretData,PullData,PullFretData

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


def constants(**kwargs):
  for (key,val) in kwargs.iteritems():
    if hasattr(Base,key):
      pass

class Container(object):
  pass

class Base(object):
  ".fret .f .ext and other meta-data (sample rate, pull speeds, )"
  # also classmethods which change experimental constants like bead radii,
  # trap stiffnesses, etc.
  def __init__(self,*args):#,fret,fext):
    self.fieldnames = ()
    for arg in args:
      self._setattr(arg)

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
  # pull[0],pull[1], etc. acts like a list with each being different pull
  # pull.molname = 's1m4'
  # pull.construct = 'SJF'
  # pull.conditions = '1B 1.0 nM'
  #
  # pull.plot(**kwargs) = sets up good defaults for title, names, etc.
  # pull.ext pull.f pull.sep pull.fret pull.donor

  def __init__(self,finfo,pulls,fret=[None]):
    if not isinstance(finfo,list):
      finfo = [finfo]
      pulls = [pulls]

    self.size = len(pulls)

    if fret != [None] and len(fret) != self.size:
      raise ExperimentError("Must have equal number of pulling curves and fret traces")

    super(Pulling,self).__init__(finfo[0])#,fret,pulls)
    self.pull = None

    self._pulls = pulls
    self._fret = fret
    self.info = finfo

  def __getitem__(self,key):
    if self._fret[key] is not None:
      return PullFretData(*(self._pulls[key]+self._fret[key]))
      #return Base(self.info[key], self._pulls[key], self._fret[key])
    else:
      return PullData(*self._pulls[key])

  def __len__(self):
    return self.size

  def __iter__(self):
    for i in range(self.size):
      yield self[i]

  @property
  def molname(self):
    return 's{}m{}'.format(self.slide,self.mol)

  @property
  def fret(self):
    return np.concatenate(map(operator.attrgetter('fret'),self._fret))

  def _concat_property(self, name):
    return map( operator.attrgetter(name), getattr(self,name))

class OpenLoop(Base):
  "camera .cam image. img"
  pass

class ForceClamp(Base):
  "force-extension .tdms camera .cam image .img"
  pass
