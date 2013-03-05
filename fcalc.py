import operator
from itertools import izip

import matplotlib.pyplot as plt
from numpy import concatenate

import image
import fileIO
import constants
from datatypes import FretData

molID = lambda t: 's{0}m{1}'.format(t.slide,t.mol)
molname = lambda t: 's{0}m{1}_{2}'.format(t.slide,t.mol,t.pull)
pN = lambda f: 'f'+str(f)+'pN'

BACKGROUND = 'background'

def info(s):
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
	verbose=True, ext=fileIO.FRET_FILE):
  "Process files given in flist as images"

  if background:
    BG = image.fromFile(background,background=True)
  else:
    BG = constants.default_background_subtract

  if isinstance(roi,str):
    roi = image.ROI.fromFile(roi)

  for fname in flist:
    if verbose: info('Opening %s...' % fname)
    img = image.fromFile(fname) - BG
    img.addROI(*roi)
    output = calculate(img)
    if verbose: info('Saving .fret data to file...')
    toFile(fileIO.change_extension(fname,ext), output, img.metadata)

def calcToFile(stack, filename, **kwargs):
  "Calculates and saves saves donor, acceptor, calculated FRET values to 3 column text file"
  fretdata = calculate(stack, **kwargs)
  toFile(filename, fretdata, stack.metadata)
  return fretdata

def calculate(stack, beta=constants.beta, gamma=constants.gamma, minsub=False):
  """Calculates FRET of a pull from an image.Stack

  calculate( image.Stack, beta = constants.beta, gamma = constants.gamma)

  RETURNS array of calculated FRET for each frame
  """
  donor = stack.donor - (minsub and min(stack.donor))
  acceptor = stack.acceptor - donor*beta
  acceptor = acceptor - (minsub and min(acceptor))
  return FretData(stack.time, donor, acceptor, acceptor/(acceptor+gamma*donor))

def toFile(filename, data, metadata, comments=''):
  return fileIO.savefret(filename, data, metadata, comments)

def fromFile(filename, **kwargs):
  return fileIO.load(filename, **kwargs)
