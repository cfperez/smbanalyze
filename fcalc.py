import numpy as np
import image
import fileIO
import constants
from fancydict import nesteddict
from datatypes import FretData

molID = lambda t: 's{0}m{1}'.format(t.slide,t.mol)
molname = lambda t: 's{0}m{1}_{2}'.format(t.slide,t.mol,t.trap)
pN = lambda f: 'f'+str(f)+'pN'

BETA = 0.0
GAMMA = 1.0

def info(s):
  print s

def warning(s):
  print s

BACKGROUND_FILENAME_FLAG = 'background'

def forster(d,R0):
    return 1/(1+(d/float(R0))**6)

def distance(e,R0):
    return R0*(1/float(e)-1)**(1./6)

def fret_model(d, R0=6.4, e0=1., e_min=0.):
    return (e0-e_min)*forster(d,R0) + e_min

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
  if isinstance(roi,str):
    roi = image.ROI.fromFile(roi)
  if len(roi) < 2:
    raise ValueError('Need 2 ROIs; given %s' % repr(roi))
  don_roi,acc_roi = roi

  BG = image.fromBackground(background)
  for fname in flist:
    try:
      if verbose: info('Opening %s...' % fname)
      img = image.fromFile(fname, roi=roi, background=BG)
      calculated_fret = fret_from_image(img, **calcOptions)
      if verbose:
        info('Saving .fret data to file...')
      calculated_fret.metadata.update(img.metadata)
      save_fret(calculated_fret, fileIO.change_extension(fname,ext))
    except IOError as e:
      warning("Error processing file {0}: {1}".format(
        fname, e.strerror))
    except image.StackError as e:
      warning("\n** Error processing file {}:\n\t{!s}**\n".format(fname, e))

def counts_from_image(img, roi, bg_pixel=0):
  return img.counts(roi) - bg_pixel*roi.size

def assert_all_listy(*args):
  for arg in args:
    assert isinstance(arg, (tuple, list))

def fret_from_image(img, bg_pixel=(0,0), bg_counts=(0,0), beta=BETA, gamma=GAMMA, rois=()):
  '''Returns (time,donor,acceptor,fret) calculated from an image with added ROIs

  img         image.Stack() or equivalent. Must have .roi dictionary with 'donor' and 'acceptor'
  bg_pixel    tuple of counts per pixel to substract from (donor,acceptor)
  bg_counts   tuple of counts to substract from (donor,acceptor)
  '''
  assert_all_listy(bg_pixel, bg_counts, rois)
  if not img.roi:
    if not rois:
      raise ValueError('Must specify ROIs in image or provide using rois=')
    else:
      img.addROI(*rois)
  pixel_donor, pixel_acc = bg_pixel
  bg_donor, bg_acc = bg_counts
  donor_roi,acceptor_roi = img.roi['donor'], img.roi['acceptor']
  don = counts_from_image(img, donor_roi, pixel_donor) - bg_donor
  acc = counts_from_image(img, acceptor_roi, pixel_acc) - bg_acc
  dcounts,acounts,fret = fret_counts(don, acc, beta, gamma)
  metadata = nesteddict.from_dict(img.metadata)
  metadata.update(
    beta=beta, gamma=gamma, bg_pixel=bg_pixel, bg_counts=bg_counts,
    roi_donor=donor_roi.toDict(), roi_acceptor=acceptor_roi.toDict()
    )
  fdata = FretData.fromFields(img.time,dcounts,acounts,fret)
  fdata.metadata = metadata
  return fdata

def fret_counts(donor, acceptor, beta=BETA, gamma=GAMMA):
  """Returns (donor,acceptor,fret) adjusted for crosstalk beta, gamma
  """
  beta,gamma = float(beta),float(gamma)
  acceptor_ = acceptor - donor*beta
  donor_ = donor
  return donor_, acceptor_, acceptor_/(acceptor_+gamma*donor_)

def calculate(donor, acceptor, beta=BETA, gamma=GAMMA):
  """Returns (donor,acceptor,fret) adjusted for crosstalk beta, gamma
  """
  print "Use fret_counts() instead"
  return fret_counts(donor,acceptor,beta,gamma)

def calc_gamma(pre, post):
    diff_acc = np.mean(pre.acceptor) - np.mean(post.acceptor)
    diff_don = np.mean(post.donor) - np.mean(pre.donor)
    return diff_acc/diff_don

def save_fret(fret, filename, comments=''):
  return fileIO.savefret(filename, fret.time,
   fret.donor, fret.acceptor, fret.fret,
    fret.metadata, comments)

def fromFile(filename, **kwargs):
  return FretData.fromFile(filename)