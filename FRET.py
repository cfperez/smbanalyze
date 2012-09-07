import Image, FileIO
import matplotlib.pyplot as plt
from FileIO import savedat
import os, glob, re

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

def fromDirectory(dir=None, verbose=False, **kwargs):

    old_dir = os.getcwd()
    if dir:
	os.chdir(dir)

    files = glob.glob( '*.img' )
    bg_files = glob.glob( '*_background.img' )

    if verbose:
	print "Found files:\n%s\n" % '\n'.join(files)

    background = ''
    BG = None
    results = {}

    for file in files:
	
	# skip background files for FRET processing
	if file in bg_files:
	    continue

	basename, ext = os.path.splitext(file)

	def find_bg_filename(basename):
	    if basename == '':
		return ''
	    bgName = basename+'_background.img'
	    if bgName in bg_files: 
		return bgName if bgName.count('_') > background.count('_') else background
	    else:
		return find_bg_filename(basename.rpartition('_')[0])

	background = find_bg_filename(basename)

	if verbose:
	    print "Processing image %s using background %s" % (file,background)

	image = Image.Stack(file)

	if background:
	    if BG is None or background != BG.filename:
		BG = Image.fromBackground(background)
	    image = image - BG

	data = calctofile( image, basename+'.fret' )

	pattern = re.compile(r'^([a-zA-Z0-9]+)_.*s(\d+)m(\d+)(?:_(\d+))?(?:_(background))?$')
	construct, slide, mol, pull, isBackground = pattern.match(basename).groups()

	temp = Experiment()
	temp.image = image
	temp.fret = data
	results['s%sm%sp%s'%(slide,mol,pull)] = temp

	if kwargs.get('plotall') is True:
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
    return results

class Experiment:
    pass
