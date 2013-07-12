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

BACKGROUND = 'background'

def info(s):
  print s

def warning(s):
  print s

def processMatch(*fglob, **kwargs):
  '''Process image files found using fmatch() of arguments

  Convenience function of processFiles(), which gets passed all the
  keyword arguments. Check processFiles() docs for options.
  '''
  fglob = fglob + (fileIO.IMAGE_FILE,)
  filelist = filter(lambda s: s.count(BACKGROUND)==0, fileIO.fmatch(*fglob))
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
      output = calculate(img, **calcOptions)
      if verbose: info('Saving .fret data to file...')
      toFile(fileIO.change_extension(fname,ext), output, img.metadata)
    except IOError as e:
      warning("Error processing file {0}: {1}".format(
        fname, e.strerror)
      )

def calcToFile(stack, filename, **kwargs):
  "Calculates and saves saves donor, acceptor, calculated FRET values to 3 column text file"
  fretdata = calculate(stack, **kwargs)
  toFile(filename, fretdata, stack.metadata)
  return fretdata

def calculate(stack, beta=constants.beta, gamma=constants.gamma, minsub=False):
  """Calculates FRET from an image.Stack

  calculate( image.Stack, beta = constants.beta, gamma = constants.gamma)

  RETURNS array of calculated FRET for each frame
  """
  donor = stack.donor - (minsub and min(stack.donor))
  acceptor = stack.acceptor - donor*beta
  acceptor = acceptor - (minsub and min(acceptor))
  return FretData.fromFields(stack.time, donor, acceptor, acceptor/(acceptor+gamma*donor))

def toFile(filename, data, metadata, comments=''):
  return fileIO.savefret(filename, data, metadata, comments)

def fromFile(filename, **kwargs):
  return FretData.fromFile(filename)
