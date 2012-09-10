import Image, FileIO
import matplotlib.pyplot as plt
from FileIO import savedat
import os, glob, re, useful

beta = 0.13
gamma = 1.16

def calc(stack, beta=beta, gamma=gamma):
    """Calculates FRET of a pull from an Image.Stack

calcFRET( Image.Stack, beta = Image.beta, gamma = Image.gamma)

defaults: taken from previous measurements

RETURNS array of calculated FRET for each frame
"""

    donor = stack.donor - min(stack.donor)
    acceptor = stack.acceptor - donor*beta
    acceptor = acceptor - min(acceptor)

    return acceptor/(acceptor+gamma*donor)

def calctofile(stack, filename, **kwargs):
    "saveFRETdata( fret, ImageStack, filename): saves donor,acceptor, FRET to 3 column text file"

    fretdata = calc(stack, **kwargs)
    savedat(filename, (stack.donor,stack.acceptor,fretdata), header='donor acceptor FRET', fmt=('%u','%u','%.5f'))
    return fretdata

def fromDirectory(*args, **kwargs):

    dir = kwargs.get('dir')
    verbose = 'verbose' in args or kwargs.get('verbose')
    plotall = 'plotall' in args or kwargs.get('plotall')
    roi_file = kwargs.get('roi_file') # 'roi*' would be the convention

    old_dir = os.getcwd()
    if dir:
	os.chdir(dir)

    roi_file = glob.glob(roi_file) if roi_file else None
    if roi_file:
	roi = roi_file.pop()
	Image.setDefaultROI( \
	    *Image.ROI.fromfile(roi, origin=kwargs.get('roi_origin','absolute')))

	if roi_file:
	    print "WARNING: Only using first ROI file found: %s" % roi
	elif verbose:
	    print "Using ROI file: %s" % roi
    elif not Image.Stack.defaultROI:
	raise RuntimeError, "No ROIs set or loaded--cannot compute counts"

    files = glob.glob( '*.img' )
    bg_files = glob.glob( '*_background.img' )

    if verbose:
	print "Found files:\n%s\n" % '\n'.join(files)

    background = ''
    slide_background = ''
    last_slide = None
    BG = None
    results = Experiments()

    pattern = re.compile(FileIO.FILENAME_SYNTAX)

    for file in files:
	
	# skip background files for FRET processing
	if file in bg_files:
	    continue

	basename, ext = os.path.splitext(file)

	construct, slide, mol, pull, time, series, isBackground = \
	    pattern.match(basename).groups()

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
	    background = bg_search
	elif slide == last_slide and bg_search.count('_') < background.count('_'):
	    background = slide_background

	image = Image.Stack(file)

	if background:
	    if BG is None or background != BG.filename:
		BG = Image.fromBackground(background)
	    image = image - BG

	if verbose:
	    print "Processing image %s using background %s" % (file,background)

	data = calctofile( image, basename+'.fret' )


	pull = int(pull) or 1

	temp = Experiment(image=image, fret=data)

	molID = 's%sm%s'%(slide,mol)
	results[construct][molID][pull] = temp

	last_slide = slide

	if plotall:
	    plt.figure()
	    plt.subplot(211)
	    plt.title(' '.join([construct, 's%sm%sp%s'%(slide,mol,pull)]) )
	    plt.plot(image.donor, label='donor')
	    plt.plot(image.acceptor,'r-', label='acceptor')
	    plt.ylabel('counts')
	    plt.legend(loc='upper left')

	    plt.subplot(212)
	    plt.ylabel('FRET')
	    plt.plot(data, 'g-',label='fret')
    
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
