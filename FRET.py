import Image, FileIO
import matplotlib.pyplot as plt
from numpy import concatenate
from FileIO import savedat,loaddat
import os, glob, re, useful, operator

beta = 0.13
gamma = 1.16

molname = lambda t: 's{}m{}_{}'.format(t.slide,t.mol,t.pull)

def plot(*data, **kwargs):
  "loc=legend location (or None for off), title, hold (False for individual plots),\
  names=labels for individual traces (must equal number of traces)\
  prefix=for all traces in legend"

  loc = kwargs.get('legend', 'best')
  title = kwargs.get('title')
  hold = kwargs.get('hold',True)

  # to prefix labels in plots
  names = kwargs.get('names',[None]*len(data))
  prefix = kwargs.get('prefix','')

  if len(names) != len(data):
	raise ValueError, "Must have same number of names as traces to plot"

  for name,trace in zip(names,data):
	hasfret = kwargs.get('fret',hasattr(trace,'fret'))

	if not hold:
	  plt.figure()

	if hasfret:
		plt.subplot(211)

	if title:
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
	  plt.legend(loc=loc)

	if hasfret:
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

  counts,bins,patches = \
	plt.hist( concatenate(map(operator.attrgetter(attr),data)), bins )

  delta = bins[1]-bins[0]
  bins = (bins-delta/2)[1:]

  return counts,bins

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
    savedat(filename, (stack.time,stack.donor,stack.acceptor,fretdata), header='time donor acceptor FRET', fmt=('%u','%u','%.5f'))
    return fretdata

def fromDirectory(*args, **kwargs):

	dir = kwargs.get('dir')
	# roi_origin => for roi file
	verbose = 'verbose' in args or kwargs.get('verbose')
	plotall = 'plotall' in args or kwargs.get('plotall')
	hold = 'singleplot' in args or kwargs.get('singleplot')
	saveplot = 'saveplot' in args or kwargs.get('saveplot')
	roi_file = kwargs.get('roi_file') # 'roi*' would be the convention
	files = kwargs.get('glob','*.img')
	user_slide = kwargs.get('slide')
	user_mol = kwargs.get('mol')
	background = kwargs.get('background','')

	old_dir = os.getcwd()
	if dir:
	  os.chdir(dir)

	try:
	  roi_file = glob.glob(roi_file) if roi_file else None
	  if roi_file:
		roi = roi_file.pop()
		Image.setDefaultROI(
		  *Image.ROI.fromFile(roi, origin=kwargs.get('roi_origin','absolute')))

	  if roi_file:
		print "WARNING: Only using first ROI file found: %s" % roi_file
	  elif verbose:
		print "Using ROI file: %s" % roi_file
	  elif not Image.Stack.defaultROI:
		raise RuntimeError, "No ROIs set or loaded--cannot compute counts"

	  files = glob.glob( '*.img' )
	  bg_files = glob.glob( '*_background.img' )

	  if verbose:
		print "Found files:\n%s\n" % '\n'.join(files)

	  slide_background = ''
	  last_slide = None
	  BG = None
	  results = Experiments()

	  for fname in files:

	  # skip background files for FRET processing
		if fname in bg_files:
		  continue

		basename, ext = os.path.splitext(fname)

		finfo = FileIO.parseFilename(fname)
		construct, context, slide, mol, pull, force, min, sec,\
		  series, isBackground = finfo

		def match_or_included(x,y):
		  if x is not None:
			if type(x) is tuple:
			  return y in x
			else:
			  return x==y
		  return True

		if not (match_or_included(user_slide,slide) or \
			match_or_included(user_mol,mol)):
		  print "skipping " + fname
		  continue

		# recursively search for background file with the most specific scope
		# using '_' convention of file naming
		def find_bg_filename(basename):
			if basename == '':
			  return ''
			bgName = basename+'_background.img'
			if bgName in bg_files: 
			  return bgName if bgName.count('_') >= background.count('_') \
				else background
			else:
			  return find_bg_filename(basename.rpartition('_')[0])

		bg_search = find_bg_filename(basename)
		if slide != last_slide or bg_search.count('_') >= background.count('_'):
		  slide_background = background
		  background = bg_search or background
		elif slide == last_slide and bg_search.count('_') < background.count('_'):
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

		temp = Experiment(image=image, fret=data)

		molID = 's%sm%s'%(slide,mol)
		results[construct][molID][pull] = temp

		last_slide = slide

		if plotall:
		  figure()
		  title = ' '.join([construct, 's%sm%sp%s'%(slide,mol,pull)])
		  plot(image, fret=data, title=title)

		  if saveplot:
			plt.savefig('%s %s s%sm%s_%s.png'%(construct,context.replace('_',' '), \
				slide,mol,pull) )
	  
	finally:
	  os.chdir(old_dir)

	return results.__lock__()

class Exp2:
    def __init__(self,*args,**kwargs):
	self._mydict = dict(*args,**kwargs)

    def __getitem__(self,name):
	if name not in self._mydict and '=' not in name:
	    self._mydict[name] = Exp2()
	return self._mydict[name]

    def __setitem__(self,name,val):
	self._mydict[name]=val

    def __setattr__(self,name,val):
	if not name.startswith('_'):
	    try:
		a,b = re.split(r'[=,]',name)
		name = name+a
		if b:
		    value=b
	    except ValueError:
		pass
	    finally:
		self.__setitem__(name,val)
	else:
	    self.__dict__[name] = val

    def __getattr__(self,name):
	if name in self._mydict:
	    return self[name]
	else:
	    try:
		return self.__dict__[name]
	    except KeyError:
		if name.startswith('_'):
		    raise AttributeError, name
		return self[name]

    def __repr__(self):
	return repr(self._mydict)


class Experiments(dict):
    def __init__(self,*args,**kwargs):
	super(self.__class__,self).__init__(*args,**kwargs)
	self._lock = False

    def __getitem__(self,name):
	if not self._lock:
	    if type(name) == str and \
		    name not in self and '=' not in name and '[' not in name:
		self[name] = Experiments()
	return dict.__getitem__(self, name)

    def __getattr__(self, name):
	if not name.startswith('_'):
	    return self.__getitem__(name)
	else:
	    dict.__getattribute__(self,name)

    def __setattr__(self, name, value):
	if not name.startswith('_'):
	    self[name] = value
	else:
	    super(self.__class__,self).__setattr__(name,value)

    def __lock__(self):
	self._lock = True
	for item in self.items():
	    if hasattr(item,'__lock__'):
		item.__lock__()
	return self

class Experiment:
    def __init__(self, **kwargs):
	self.__dict__.update(kwargs)

    def __setattr__(self, name, value):
	raise AttributeError, "Instance of class %s is read-only" % self.__class__

    def __repr__(self):
	return 'Experiment: ' + str(self.__dict__)
