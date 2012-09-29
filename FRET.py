import os
import glob
import re
import operator

import matplotlib.pyplot as plt
from numpy import concatenate,iterable

import useful
import Image
import FileIO

beta = 0.13
gamma = 1.16

molID = lambda t: 's{}m{}'.format(t.slide,t.mol)
molname = lambda t: 's{}m{}_{}'.format(t.slide,t.mol,t.pull)
pN = lambda f: 'f'+str(f)+'pN'

def plot(*data, **kwargs):
  "loc=legend location (or None for off), title, hold (False for individual plots),\
  names=labels for individual traces (must equal number of traces)\
  prefix=for all traces in legend"

  loc = kwargs.get('legend', 'best')
  hold = kwargs.get('hold',True)

  if len(data) == 1 and iterable(data[0]):
	data = data[0]

  # to prefix labels in plots
  names = kwargs.get('names',[None]*len(data))
  prefix = kwargs.get('prefix','')
  title = kwargs.get('title','')

  if len(names) != len(data):
	raise ValueError, "Must have same number of names as traces to plot"

  for name,trace in zip(names,data):
	hasfret = kwargs.get('fret',hasattr(trace,'fret'))

	if not hold:
	  plt.figure()

	if hasfret is not False:
		plt.subplot(211)

	plt.title(title)

	name = name or getattr(trace,'molname','')

	x_axis = None
	if hasattr(trace,'time'):
	  x_axis = trace.time
	  plt.xlabel('Seconds')
	else:
	  x_axis = range(1,len(trace.donor)+1)
	  plt.xlabel('Frames')

	plt.plot(x_axis, trace.donor, label=' '.join([prefix,name,'donor']))
	plt.plot(x_axis, trace.acceptor,label=' '.join([prefix,name,'acceptor']))
	plt.ylabel('Counts')
	if loc:
	  plt.legend(loc=loc,ncol=2,prop={'size':'small'})

	if hasfret is not False:
		plt.subplot(212)
		plt.ylabel('FRET')
		plt.xlabel('Frames')
		if hasfret is True:
			plt.plot( trace.fret, label=' '.join([prefix,name]))
		else:
			plt.plot(hasfret, label=' '.join([prefix,name]))

def hist(*data, **kwargs):
  bins = kwargs.get('bins',50)
  attr = kwargs.get('plot','fret')

  #if len(data) == 1:
  #  data = list(data[0])

  if kwargs.get('join',True):
    func = concatenate
  else:
    func = lambda x: x

  counts,bins,patches = plt.hist(
        func(map(operator.attrgetter(attr),data)), bins)

  delta = bins[1]-bins[0]
  bins = (bins-delta/2)[1:]

  return bins,counts

def calc(stack, beta=beta, gamma=gamma):
    """Calculates FRET of a pull from an Image.Stack

calcFRET( Image.Stack, beta = Image.beta, gamma = Image.gamma)

RETURNS array of calculated FRET for each frame
"""

    donor = stack.donor - min(stack.donor)
    acceptor = stack.acceptor - donor*beta
    acceptor = acceptor - min(acceptor)

    return acceptor/(acceptor+gamma*donor)

def calcToFile(stack, filename, **kwargs):
    "saveFRETdata( fret, ImageStack, filename): saves donor,acceptor, FRET to 3 column text file"

    fretdata = calc(stack, **kwargs)
    FileIO.savedat(filename, (stack.time,stack.donor,stack.acceptor,fretdata), header='time donor acceptor FRET', fmt=('%.3f','%u','%u','%.5f'))
    return fretdata

def calcDirectory(file_glob='',*args, **kwargs):

	dir = kwargs.get('dir')
	# roi_origin => for roi file
	verbose = 'verbose' in args or kwargs.get('verbose')
	plotall = 'plotall' in args or kwargs.get('plotall')
	hold = 'singleplot' in args or kwargs.get('singleplot')
	saveplot = 'saveplot' in args or kwargs.get('saveplot')
	roi_file = kwargs.get('roi_file','roi.txt') # 'roi*' would be the convention
	user_slide = kwargs.get('slide')
	user_mol = kwargs.get('mol')
	background = kwargs.get('background','')

	old_dir = os.getcwd()
	if dir:
	  os.chdir(dir)
    
	try:
	  if os.path.isfile(roi_file):
		Image.setDefaultROI(
		  *Image.ROI.fromFile(roi_file, origin=kwargs.get('roi_origin','absolute')))
		if verbose:
		  print "Using ROI file: %s" % roi_file

	  if not Image.Stack.defaultROI:
		raise RuntimeError, "No ROIs set or loaded--cannot compute counts"

	  files = glob.glob('*'+file_glob+'*.img')
	  bg_files = glob.glob( '*_background.img' )

	  if verbose:
		print "Found files:\n%s\n" % '\n'.join(files)

	  slide_background = '' # last background used on this slide
	  last_slide = None     # slide number of last file
	  BG = None             # current background image
	  results = useful.dotdict() # output

	  for fname in files:

        # Skip background files for FRET processing
		if fname in bg_files:
		  continue

		basename, ext = os.path.splitext(fname)

        # Get molecule information from filename
		finfo = FileIO.parseFilename(fname)
		(construct, context, slide, mol, pull, force, min, sec,
            series, isBackground) = finfo

        # Only process requested slide and molecules
		if (not match_or_included(user_slide,slide) or 
		  not match_or_included(user_mol,mol)):
		  if verbose:
			print "skipping " + fname
		  continue

		def find_bg_filename(basename):
			if basename == '':
			  return ''
			bgName = basename+'_background.img'
			if bgName in bg_files: 
			  return bgName if bgName.count('_') >= background.count('_') \
				else background
			else:
			  return find_bg_filename(basename.rpartition('_')[0])

        # Look for background file using name heirarchy
		bg_search = find_bg_filename(basename)

        # Use the new background if a) on a new slide or 
		# b) it is more appropriate
		if slide != last_slide or bg_search.count('_') >= background.count('_'):
		  slide_background = background
          # keep the last background if a new one is not found
		  background = bg_search or background
        # Use the last when when, e.g. s1m1_background and then s1m2
		elif slide == last_slide and \
		   bg_search.count('_') < background.count('_'):
		  background = slide_background

		image = Image.Stack(fname)

		if background:
		  if BG is None or background != BG.filename:
			BG = Image.fromBackground(background)
		  image = image - BG

		if verbose:
		  print "Processing image file '%s' using background '%s'" \
			% (fname,background)

		data = calcToFile( image, basename+'.fret' )

		temp = useful.dotdict(image=image, fret=data, 
		  donor=image.donor, acceptor=image.acceptor,molname=molname(finfo))
		if force is not None:
		  results[construct][molID(finfo)][pN(force)][series] = temp
		else:
		  results[construct][molID(finfo)][pull] = temp

		last_slide = slide

		if plotall:
		  plt.figure()
		  plot(image, fret=data, title=' '.join([construct, molname(finfo)]))

		  if saveplot:
			plt.savefig('%s %s s%sm%s_%s.png' % (construct,
			  context.replace('_',' '), slide,mol,pull) )
	  
	finally:
	  os.chdir(old_dir)

	results._lock()
	return results

# recursively search for background file with the most specific scope
# using '_' convention of file naming
def match_or_included(x,y):
  if x is not None:
    if isinstance(x,tuple):
      return y in x
    else:
      return x==y
  return True
