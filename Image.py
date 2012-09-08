from __future__ import with_statement
import numpy as np
import matplotlib.pyplot as plt
import os.path, copy, useful, FileIO

class StackError(Exception):
    pass

class ROIError(Exception):
    pass

def setDefaultROI(*args):
    Stack.setDefaultROI(*args)
    return Stack

def fromBackground(filename, filter='median'):
    bg = Stack(filename)
    width, height = bg.width, bg.height
    if filter == 'median':
	bg._img = np.median( bg._img, axis=0, overwrite_input=True ).reshape((1,height,width))
    elif filter == 'mean':
	bg._img = np.mean( bg._img, axis=0 ).reshape((1,height,width))
    elif filter == 'min':
	bg._img = np.min( bg._img, axis=0 ).reshape((1,height,width))
    else:
	raise ValueError, "Filter type can only be median, mean, or min"
	
    return bg
    
################################################################################
##
## Class ROI => Region of Interest for extracting data from images
##
## 
################################################################################
class ROI:
    """Stores static ROI data

name => how you will call the ROI after attached to an Image.Stack
origin => is the origin 'relative' to the subimage or 'absolute' to the CCD origin

Usage: 

1) Specify corners

    ROI( left, right, bottom, top, name='', origin='relative')

2) Specify bottom left and top right corner

    ROI( bottomleft, topright, name='', origin='relative' )

3) Copy from another ROI

    ROI( roi, name='', origin='relative' )
    """

    def __init__(self, left=None, right=None, bottom=None, top=None, name='', origin='relative'):

	#################################################
	## If the first argument is another ROI, make a copy
	#################################################
	try:
	    roi=left
	    if name:
		roi.name = name
	    if not hasattr(roi, 'origin'):
		roi.origin = origin

	    self.__myinit__( roi.left, roi.right, roi.bottom, roi.top, roi.name, roi.origin )
	except AttributeError:
	    
	    #################################################
	    ## If the first argument is a tuple, assume
	    ## bottomleft and topright (bounding box) usage
	    #################################################
	    try:
		self.__myinit__(left[0], right[0], left[1], right[1], name, origin)
	    except TypeError, IndexError:

		#################################################
		## Otherwise, use all the corners
		#################################################
		self.__myinit__(left,right,bottom,top,name,origin)

    def __myinit__(self, left, right, bottom, top, name, origin):
	if left is not None:
	    if right < left or left < 0 or right < 0:
		raise ROIError, "ROI must satisfy condition 0 <= left < right: left=%d right=%d" % (left,right)
	    if bottom > top or bottom < 0 or top < 0:
		raise ROIError, "ROI must satisfy condition 0 <= bottom < top: bottom=%d top=%d" % (bottom,top)

	self.left = left
	self.right = right
	self.top = top
	self.bottom = bottom
	self.name = name
	self.origin = origin

    @classmethod
    def fromfile(cls, filename):
	toInt = lambda s: int(useful.toNum(s))
	settings = FileIO.loadsettings(filename, cast=toInt)

	self = []
	for name, roi in settings.items():
	    self.append( cls(roi, name=name, origin='absolute') )

	return tuple(self)

    def toRelative(self,first,*args):
	pass
	    
    @property
    def width(self):
	return self.right-self.left+1

    @property
    def height(self):
	return self.top-self.bottom+1

    def __repr__(self):
	if None in (self.left, self.right, self.bottom, self.top):
	    name = self.name or 'Undefined'
	    return "ROI '%s' = uninitialized" % self.name
	else:
	    return "ROI '%s' = L: %d, R: %d, B: %d, T: %d" % \
		(self.name,self.left,self.right,self.bottom,self.top) 

################################################################################
##
## Class Stack
##
## 
################################################################################
class Stack:
    "A class to load and manipulate images generated during single molecule experiments"

    defaultROI = {}
    defaultDonorROI = 'donor'
    defaultAcceptorROI = 'acceptor'

    def __init__(self, filename, filetype="img", camFile='', deepcopy=False):
	if filetype != 'img':
	    raise ValueError, "Only filetype 'img' is supported!"

	#################################################
	## Load data from file called filename if string
	#################################################
	if isinstance(filename, str):
	    self.filename = filename
	    self._img = FileIO.loadimg(filename)
	    self._roi = {}
	    self._donorROIName = Stack.defaultDonorROI
	    self._acceptorROIName = Stack.defaultAcceptorROI

	    self.addROI(*self.__class__.defaultROI.values())

	    camFile = camFile or (os.path.splitext(filename)[0] + '.cam')
	    
	    settings = FileIO.loadcam(camFile)
	    self._settings = []
	    for (setting, value) in settings.iteritems():
		if not hasattr(self,setting):
		    setattr(self, setting, value)
		    self._settings += [setting]

	    # check cam and img file correspondence
	    if self._img.shape != (settings['frames'],settings['height'],settings['width']):
		raise StackError, ".img file and .cam file dimensions do not agree"

	    self._origin = (self.roileft,self.roibottom)

	#################################################
	## Make a copy of filename if actually another Stack
	#################################################
	elif isinstance(filename, Stack):
	    if deepcopy:
		self._img=filename._img.copy()
	    else:
		self._img=filename._img

	    self._roi=filename._roi
	    self._donorROIName = filename._donorROIName
	    self._acceptorROIName = filename._acceptorROIName
	    self._origin = filename._origin
	    self._settings = filename._settings
	    for setting in self._settings:
		if not hasattr(self,setting.lower()):
		    setattr(self, setting, getattr(filename,setting))

	else:
	    raise StackError, "Invalid constructor call using %s" % str(filename)

    @property
    def frames(self):
	return self._img.shape[0]

    @property
    def height(self):
	return self._img.shape[1]

    @property
    def width(self):
	return self._img.shape[2]

    @property
    def donor(self):
	if not self._roi.has_key(self._donorROIName):
	    raise StackError, "ROI called %s hasn't been defined yet" % self._donorROIName

	return self.counts(self._roi[self._donorROIName])

    @property
    def acceptor(self):
	if not self._roi.has_key(self._acceptorROIName):
	    raise StackError, "ROI called %s hasn't been defined yet" % self._acceptorROIName

	return self.counts(self._roi[self._acceptorROIName])
	

    @classmethod
    def setDefaultROI(cls, *args):
	for _roi in args:
	    cls.defaultROI[_roi.name]=_roi

    def addROI(self, *ROIs):
	for roi in ROIs:
	    try:
		key = roi.name

		if roi.origin == 'absolute':
		    # recast to relative origin
		    roi = ROI( roi.left-self.roileft, roi.right-self.roileft, \
			roi.bottom-self.roibottom, roi.top-self.roibottom, name=key )
		else:
		    # make a new object
		    roi = ROI(roi)

		if roi.right > self.width or roi.top > self.height:
		    raise ROIError

	    except AttributeError:
		raise TypeError, "Must use objects with ROI interface"
	    except ROIError:
		raise ROIError, "ROI out of bounds! %s" % repr(roi)

	    self._roi[key] = roi

	return self

    def setDonorROI(self, roi_name):
	if not self._roi.has_key(roi_name):
	    raise KeyError, "Image.Stack does not have an ROI named %s" % roi_name

	self._donorROIName = roi_name
	return self

    def setAcceptorROI(self, roi_name):
	if not self._roi.has_key(roi_name):
	    raise KeyError, "Image.Stack does not have an ROI named %s" % roi_name

	self._acceptorROIName = roi_name
	return self
	
    def counts(self, roi=None):
	if roi:
	    return self[:,roi.bottom:roi.top,roi.left:roi.right].counts()
	    #return np.sum( np.sum(self._img[:,roi.bottom:roi.top,roi.left:roi.right],axis=1), axis=1 )
	else:
	    return np.sum( np.sum(self._img,axis=1), axis=1 )

    def __getitem__(self,key):
	if type(key) == str and self._roi.has_key(key):
	    return self.counts(self._roi[key])
	elif type(key) == int:
	    return Frame(self._img[key])
	else:
	    temp = Stack(self)
	    temp._img = temp._img[key]
	    if isinstance(temp._img.shape,tuple) and len(temp._img.shape) > 2:
		temp.frames = temp._img.shape[0]
	    else:
		temp.frames = 1

	    return temp
	raise IndexError, "Invalid index: %s" % str(key)

    #def __slice__(self, slice):
	#temp = Stack(self)
	#temp._img = temp._img[slice]
	#return temp

    def append(self, stack):
	temp = copy.copy(self)
	temp._img = np.append( temp._img, stack._img, axis=0 )
	return temp

    def __sub__(self, stack):
	return self.__add__(stack.__neg__())

    def __add__(self, stack):
	temp = Stack(self)
	temp._img = temp._img + stack._img
	return temp

    def __neg__(self):
	temp = Stack(self)
	temp._img = -temp._img
	return temp

    def __repr__(self):
	return "Stack %dx%dx%d" % (self.frames, self.height, self.width)

class Frame:
    
    def __init__(self, imgarray):
	self._img = imgarray

    def __getitem__(self,key):
	#if type(key) == tuple:
	#    return self._img[key[1],key[0]]
	if isinstance(key, ROI):
	    return self._img[ key.bottom:key.top, key.left:key.right ]
	else:
	    return self._img[:,key]

    @property
    def height(self):
	return self._img.shape[0]

    @property
    def width(self):
	return self._img.shape[1]

    def counts( roi ):
	return np.sum( self[roi] )

    def show(self):
	plot = self._plot_obj = plt.imshow(self._img)
	plot.get_axes().invert_yaxis()
	plt.draw()
	return self._plot_obj

    def __neg__(self):
	temp = copy.copy(self)
	temp._img = -temp._img
	return temp

    def __repr__(self):
	return "Image %dx%d" % self._img.shape

