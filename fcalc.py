import operator
from itertools import izip

import matplotlib.pyplot as plt
from numpy import concatenate

import image
import fileIO
import constants
from datatypes import FretData

molID = lambda t: 's{0}m{1}'.format(t.slide,t.mol)
molname = lambda t: 's{0}m{1}_{2}'.format(t.slide,t.mol,t.trap)
pN = lambda f: 'f'+str(f)+'pN'


def info(s):
  print s

def warning(s):
  print s

BACKGROUND_FILENAME_FLAG = 'background'

def forster(x,R):
    return 1/(1+(float(x)/R)**6)

def distance(e,R):
    return R*(1/float(e)-1)**(1./6)

def processMatch(*fglob, **kwargs):
  '''Process image files found using fmatch() of arguments

  Convenience function of processFiles(), which gets passed all the
  keyword arguments. Check processFiles() docs for options.
  '''
  fglob = fglob + (fileIO.IMAGE_FILE,)
  filelist = filter(lambda s: s.count(BACKGROUND_FILENAME_FLAG)==0, fileIO.flist(*fglob))
  processFiles(filelist, **kwargs)

def processFiles(flist, roi='roi.txt', background=None, 
	verbose=True, ext=fileIO.FRET_FILE, **calcOptions):
  "Calculate donor, acceptor, and FRET values for files in flist argument."

  if isinstance(background, str):
    BG = image.fromFile(background, background=True)
  elif isinstance(background, int):
    BG = background
  else:
    BG = constants.default_background_subtract

  if isinstance(roi,str):
    roi = image.ROI.fromFile(roi)

  for fname in flist:
    try:
      if verbose: info('Opening %s...' % fname)
      img = image.fromFile(fname) - BG
      img.addROI(*roi)
      donor,acc,fret = calculate(img.donor, img.acceptor, **calcOptions)
      if verbose: info('Saving .fret data to file...')
      toFile(fileIO.change_extension(fname,ext), 
        img.time, donor, acc, fret, 
        img.metadata)
    except IOError as e:
      warning("Error processing file {0}: {1}".format(
        fname, e.strerror))
    except image.StackError as e:
      warning("\n** Error processing file {}:\n\t{!s}**\n".format(fname, e))

def counts_from_image(img, roi, bg_pixel):
  return img.counts(roi) - bg_pixel*roi.size

def calculate(donor, acceptor, beta=0.0, gamma=1.0):
  """Returns (donor,acceptor,fret) adjusted for crosstalk beta, gamma
  """
  beta,gamma = float(beta),float(gamma)
  acceptor_ = acceptor - donor*beta
  donor_ = donor + (1+beta/gamma)
  return donor_, acceptor_, acceptor_/(acceptor_+gamma*donor_)

def subtract_min_background(exp_list):
  ''' SIDE EFFECT: subtracts observed background (minimum+1) from FretData in exp_list
  '''
  if exp_list.has('fret') != exp_list:
    raise ValueError('Argument exp_list must contain only experiments with FRET')
  fret = exp_list.get('fret')
  min_ = lambda attr: min(map(lambda x: min(getattr(x, attr)), fret))
  min_d, min_a = min_('donor'), min_('acceptor')
  for p in fret:
    p.donor -= min_d
    p.acceptor -= min_a

def calculate_fret_bg_subtracted(exp_list, beta, gamma):
  ''' SIDE EFFECT: changes FretData passed in to be background subtracted and beta/gamma adjusted
  '''
  assert operator.isNumberType(beta)
  assert operator.isNumberType(gamma)
  subtract_min_background(exp_list)
  return [calculate(p.fret, beta=beta, gamma=gamma) for p in exp_list]

def minimum_counts(fretdata):
  return min(fretdata.donor), min(fretdata.acceptor)

def toFile(filename, time, donor, acceptor, fret, metadata, comments=''):
  return fileIO.savefret(filename, time, donor, acceptor, fret, metadata, comments)

def fromFile(filename, **kwargs):
  return FretData.fromFile(filename)
