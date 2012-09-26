import os
import collections
import operator

import numpy as np

from useful import dotdict
import FileIO
import FRET

parameters = {'T':293.2}

kT = lambda T: 0.0138965*T

FretData = collections.namedtuple('FretData',('time','donor','acceptor','fret'))
PullData = collections.namedtuple('PullData', ('ext','f','sep'))

class ExperimentError(Exception):
  pass

def load(filename, **kwargs):
  "Create a new Experiment subclass based on filetype and data"
  verbose = kwargs.get('verbose')
  recalc = kwargs.get('recalc')

  basename,ext = os.path.splitext(filename)
  if ext == '.str':
	if verbose:
	  print "Processing %s as pulling experiment..." % basename

	fretfile = FileIO.add_fret_ext(basename)
	if recalc:
	  imgfile = FileIO.add_img_ext(basename)
	  if verbose:
		print "Recalculating fret from %s" % imgfile
	  image = Image.Stack(imgfile)
	  fret = FRET.calcToFile(image,fretfile) 
	else:
	  if verbose:
		print "Loading fret from %s" % fretfile
	  fret = FileIO.loadFRET(fretfile)

	pullfile=FileIO.add_pull_ext(basename)
	if verbose:
	  print "Loading pulling data from %s" % pullfile
	pulldata = FileIO.loadPull(pullfile)
	
	finfo = FileIO.parseFilename(filename)

	return Pulling(finfo,fret,pulldata)


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
	  self._loadattr(arg)

  def _loadattr(self,named):
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

  def __init__(self,finfo,pulls,fret=None):
	if type(finfo) is not list:
	  finfo = [finfo]
	  pulls = [pulls]

	super(Pulling,self).__init__(finfo[0])#,fret,pulls)
	self.pull = None

	self.size = len(pulls)

	if fret is not None and len(fret) != self.size:
	  raise ExperimentError, "Must have equal number of pulling curves and fret traces"

	self._pulls = pulls
	self._fret = fret
	self.info = finfo

	if fret:
	  self.datatype = \
		collections.namedtuple('PullingFretData',PullData._fields+FretData._fields)
	else:
	  self.datatype = PullData

  def __getitem__(self,key):
	if self._fret:
	  #return self.datatype(*(self.pulls[key]+self.fret[key]))
	  return Base(self.info[key], self._pulls[key], self._fret[key])
	else:
	  return self.datatype(*self.pulls[key])

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
