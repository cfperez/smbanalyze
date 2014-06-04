import os.path as opath
from operator import and_
import logging
import cPickle as pickle
import abc
from functools import total_ordering
from itertools import ifilter
from smbanalyze.date import today, to_date
from smbanalyze.sampling import downsample as dsample

from numpy import mean

import fileIO 
import constants
from datatypes import TrapData, FretData, hasTrapData, hasFretData
from fancydict import nesteddict

logger = logging.getLogger(__name__)
logger.setLevel(constants.logLevel)
logger.addHandler(constants.logHandler)

class ExperimentError(Exception):
  pass

@total_ordering
class Base(object):
  ".fret .f .ext and other meta-data (sample rate, trap speeds, )"
  __metaclass__ = abc.ABCMeta

  INFO_FIELDS = ()

  def __init__(self, trap, fret, metadata):
    if trap and not hasTrapData(trap):
      raise ValueError(
          "__init__ argument 'trap' <{}> does not have trap data".format(trap))
    if fret and not hasFretData(fret):
      raise ValueError(
          "__init__ argument 'fret' <{}> does not have fret data".format(fret))
    self.trap = trap
    self._fret = fret
    self.metadata = nesteddict.from_dict(metadata)
    self.metadata['trap'] = getattr(trap, 'metadata', nesteddict())
    self.metadata['fret'] = getattr(fret, 'metadata', nesteddict())

  @abc.abstractmethod
  def filenameMatchesType(cls, filename):
    pass

  @property
  def filename(self):
    return self.metadata.get('filename', '')

  @property
  def info(self):
    datesort = self.metadata.get('datetime', self.metadata.get('date', None))
    return (datesort,) + fileIO.split_fname(self.filename)

  @property
  def fret(self):
    return self._fret

  @fret.setter
  def fret(self, data):
    meta = data.metadata.copy()
    if self.fret:
      for k,v in self.fret.metadata.iteritems():
        meta.setdefault(k, v)
    self._fret = FretData(data.data, meta)
    self.metadata['fret.'] = self._fret.metadata

  @fret.deleter
  def fret(self):
    self._fret = None
    self.metadata['fret'] = {}

  def __lt__(self, other):
    this, that = self.info, other.info
    return this < that

  def __eq__(self, other):
    return self.info == other.info

  def __getitem__(self, key):
    if isinstance(key, str):
      return self.metadata[key]
    else:
      return self.get(key)

  def get(self, key):
    return self.trap and self.trap[key], self.fret and self.fret[key]

  def where(self, *conditions):
    return self[reduce(and_, conditions)]

  def __setitem__(self, key, val):
    self.metadata[key] = val

  @property
  def figure(self):
    return self._figure

  @classmethod
  def fromFile(cls, strfile, fretfile, metadata):
    assert strfile or fretfile
    assert isinstance(strfile, str)
    assert isinstance(fretfile, (str,type(None)))
    metadata = nesteddict.from_dict(metadata)

    trap = TrapData.fromFile(strfile)
    fret = FretData.fromFile(fretfile) if fretfile else None
    trap.metadata.setdefault('date', today())
    trap_datetime = trap.metadata['date']
    metadata.setdefault('date', to_date(trap_datetime))
    metadata.setdefault('datetime', trap_datetime)
    metadata['filename'] = fileIO.splitext(strfile)[0]
    newCls = cls(trap, fret, metadata)
    assert isinstance(newCls, cls)
    assert getattr(newCls, 'filename', None) is not None
    return newCls

  def save(self, filename=None, path=''):
    if filename is None and self.filename is None:
      raise ValueError('Please choose a filename')
    filename = filename or self.filename
    if not filename.endswith('.exp'):
      filename = filename + '.exp'
    full_path = opath.join(path, filename)
    with open(full_path, 'wb') as fh:
      pickle.dump(self, fh, protocol=2)

  def __repr__(self):
    if hasattr(self, 'filename'):
      return "<Experiment.%s from '%s'>" % (self.__class__.__name__, self.filename)
    else:
      return super(Base,self).__repr__()

  def __str__(self):
    return self.__repr__()


class Pulling(Base):
  "stretching curves of single molecule .str camera .cam image .img"

  forceOffsetRange = (740,770)
  extensionOffsetRange = (13,16)
  
  EXTENSION = fileIO.PULL_FILE
  INFO_FIELDS = ('date', 'construct', 'conditions', 
    'slide', 'mol', 'pull', 'reverse', 'collapsed' )

  def __init__(self, trap, fret=None, metadata={}):
    super(Pulling, self).__init__(trap, fret, metadata)
    self._ext_offset = 0

    for k in ifilter(lambda s: s.startswith('fret.'), trap.metadata.copy()):
      self.metadata['fret'][k[5:]] = trap.metadata.pop(k)

    sampling_time = self['trap'].get('sampling_time', 
      constants.default_pulling_sampling_time)
    step_size = self['trap'].get('step_size', None)
    if sampling_time and step_size:
      self['trap'].setdefault('pulling_rate', step_size/sampling_time)
    if fret:
      trap_rate = sampling_time*1000
      fret_rate = fret.metadata.pop('exposurems', None)
        # constants.default_fret_exposure_time_ms)
      if fret_rate:
        if not (fret_rate / trap_rate).is_integer():
          print 'Trap and FRET collection rates may not be even multiples! {} (trap) vs {} (fret)'.format(
                sampling_time, fret_rate)
        self.metadata['fret.exposure_time'] = fret_rate/1000.
        self.metadata['fret.exposurems'] = fret_rate
        self.metadata['sampling_ratio'] = int(fret_rate / sampling_time /1000.)

  def downsample(self, sampling_ratio=None):
    sampling_ratio = sampling_ratio or self['sampling_ratio']
    trap_ = dsample(self.trap.copy(), sampling_ratio)
    return Pulling(trap_, self.fret.copy(), self.metadata)

  @classmethod
  def fromFile(cls, strfile, fretfile=None, metadata={}):
    'Load stretching data and corresponding fret data from files'
    basename, ext = fileIO.splitext(strfile)
    strfile, fretfileFromBase = (fileIO.add_pull_ext(basename),
                                fileIO.add_fret_ext(basename))
    fretfile = fretfile or fretfileFromBase
    if not opath.exists(fretfile):
      fretfile = None
    return super(Pulling, cls).fromFile(strfile, fretfile, metadata)

  FILENAME_SYNTAX = ('construct', 'conditions', 'slide', 'mol', 'pull')

  @classmethod
  def filenameMatchesType(cls, filename):
    assert isinstance(filename, str)
    finfo = fileIO.parseFilename(filename)
    if not finfo:
      return False
    for attr in finfo._fields:
      val = getattr(finfo, attr)
      if attr not in cls.FILENAME_SYNTAX and val:
        return False
    return True

  @property
  def isReverse(self):
    '''Return True iff experiment is a reverse pull'''
    return self.metadata.get('reverse', False)

  def where(self, *conditions):
    # self[reduce(and_, conditions)]
    trap,fret = super(Pulling,self).where(*conditions)
    return Pulling(trap,fret,self.metadata)

  def forceOffset(self, x_range=None):
    x_range = x_range or Pulling.forceOffsetRange
    data = self.trap.select(x=x_range)
    if len(data) == 0:
      raise ExperimentError(
        '{0}: No data exists in range {1} - {2}'.format(
          str(self), *x_range))
    return mean(data.f)

  def adjustForceOffset(self, baseline=0.0, offset=None, offset_range=None):
    offset = offset or -self.forceOffset(offset_range)
    offset += baseline
    self.trap.f += offset
    self.trap.ext -= offset/self.trap.meanStiffness()
    return offset

  def extensionOffset(self, f_range=None, x_range=None):
    'Returns average extension of FEC between given forces'
    f_range = f_range or Pulling.extensionOffsetRange
    data = self.trap.select(x=x_range, f=f_range)
    if len(data) == 0:
      raise ExperimentError(
        '{0}: No data exists in f_range={1} and x_range={2}'.format(
          self, f_range, x_range)
      )
    return mean(data.ext)

  def adjustExtensionOffset(self, baseline, f_range=None, x_range=None):
    'Adjust extension to hit baseline. If f_range is not given, it is calculated from Pulling.extensionOffsetRange'
    offset = self.extensionOffset(x_range=x_range, f_range=f_range) - baseline
    self.trap.ext -= offset
    self._ext_offset = offset
    return offset

  # show_fret_default = False

  # def plot(self, style=None, **kwargs):
  #   show_fret = kwargs.setdefault('show_fret', Pulling.show_fret_default)
  #   FEC = kwargs.setdefault('FEC', True)
  #   if FEC:
  #     style = style or {'trap': ':'}
  #   kwargs.setdefault('legend', None)
  #   kwargs.setdefault('title', self.filename or '')
  #   kwargs.setdefault('label', self.filename or '')
  #   loc_x = min(self.trap.ext)+10
  #   location = list(kwargs.pop('annotate', (loc_x, 15)))
    # if self.fret:
  #     self.figure.plot(self.fret, self.trap, style=style, **kwargs)
  #   else:
  #     self.figure.plot(self.trap, style=style, **kwargs)
  #   # tk if reverse pull, reverse the x-axis of the FEC, which is last thing plotted
  #   self.figure.xlim(reverse=self.isReverse)
  #   return self.figure

  def savefig(self, filename=None, path=''):
    if not self.figure.exists:
      raise ExperimentError('No figure available for experiment {0}'.format(str(self)))
    else: 
      filename = filename or self.filename
      self.figure.toFile(filename)

  @classmethod
  def load(cls, filename):
    basename, extension = opath.splitext(filename)
    filename = basename + (extension or '.exp')
    return pickle.load(open(filename,'rb'))
