from __future__ import with_statement
from operator import methodcaller
import os.path
import copy

import numpy as np
import matplotlib.pyplot as plt

import FileIO
import FRET
import useful

################################################################################
## EXCEPTIONS
################################################################################
class StackError(Exception):
  pass

class ROIError(Exception):
  pass

################################################################################
## MODULE LEVEL FUNCTIONS
################################################################################
def setDefaultROI(*args):
  Stack.setDefaultROI(*args)

def fromFile(filename, **kwargs):
"fromFile(filename, [background='' or True]): Load Image from filename with optional background file to subtract"

  bg = kwargs.pop('background', False)
  try:
    img = Stack(filename, **kwargs)
  except IOError:
    raise IOError("File %s can't be found/loaded" % filename)

  if bg is True:
    NA,bg = FRET.matchImgFilestoBackground([filename])[0]
    if not bg:
      print "Warning! Did not find a background file for %s" % filename
  if bg:
    bg = fromBackground(bg, **kwargs)
  else:
    return img

  return img - bg

def fromBackground(filename, filter='median'):
  try:
    bg = Stack(filename)
  except IOError:
    raise IOError("Can't load background file")

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
  """Region of Interest

ROI( (bottom,left), (top,right), name='', origin='relative' )

name => how you will call the ROI after attached to an Image.Stack
origin => is the origin 'relative' to the subimage or 'absolute' to the CCD origin

>>> donorROI = ROI( (0,10), (0,10), name='donor' )

Other constructors from class methods:

  From a file
  >>> ROI.fromFile( 'roi.txt', [origin='absolute'] )

  Specify corners
  >>>  ROI.fromCorners( left, right, bottom, top, name='', origin='relative')

  Copy from another ROI
  >>>  ROI.copy( old_roi, name='', origin='relative' )

And saving:

  >>> ROI.save('roi.txt', ROI1, ROI2)

  or individually. Mode='a' appends to existing file, 'w' overwrites.
  >>> donorROI.toFile('filename.txt', [mode='w'])

Convert between 'absolute' and 'relative' coordinates:

  >>> donorROI.toRelative( (0,0) )
  >>> donorROI.toAbsolute( (0,0) )
  """

  def __init__(self, bl, tr, name='', origin='relative'):

    left,bottom = bl
    right,top = tr

    if right < left or left < 0 or right < 0:
      raise ROIError, "ROI must satisfy condition 0 <= left < right: left=%d right=%d" \
        % (left,right)

    if bottom > top or bottom < 0 or top < 0:
      raise ROIError, "ROI must satisfy condition 0 <= bottom < top: bottom=%d top=%d" \
        % (bottom,top)

    self.left = left
    self.right = right
    self.top = top
    self.bottom = bottom
    self.name = name
    self.origin = origin

    self.lines = []

  @classmethod
  def fromFile(cls, filename, origin='absolute'):
    "Load ROI(s) from a LabView config file and return tuple of Image.ROI objects"
    settings = FileIO.loadsettings(filename, cast=useful.toInt)

    temp = []
    for name, roi in settings.items():
      roi.filename = filename
      roi.name = name
      roi.origin = origin
      temp.append( cls.copy(roi) )

    return temp(self)

  @classmethod
  def copy(cls, roi, name='', origin=''):
    try:
      name = name or roi.name
      origin = origin or roi.origin or 'relative'
      return cls( (roi.left,roi.bottom), (roi.right,roi.top), name, origin )
    except AttributeError:
      raise ROIError, "Must use ROI-type object"

  @classmethod
  def fromCorners(cls, left, right, bottom, top, **kwargs):
      return cls( (left,bottom), (right,top), **kwargs )

  @classmethod
  def save(cls, filename, *ROIs):
    "save(filename, *ROIs): Save ROI(s) to a LabView config file format"
    mode = 'w'
    for roi in ROIs:
      roi.toFile(filename, mode)
      mode = 'a'

  Colors = ('r','y')
  _lastColor = 0

  def draw(self,color=None):
    if not color:
      i = ROI._lastColor
      color=ROI.Colors[i]
      ROI._lastColor = (i+1) % (len(ROI.Colors))

    self.clear()

    self.lines = [
      plt.axvline(self.left,color=color),
      plt.axvline(self.right,color=color),
      plt.axhline(self.bottom,color=color),
      plt.axhline(self.top,color=color)
      ]

  def clear(self):
    if self.lines:
      try:
        map( methodcaller('remove'), self.lines )
        plt.gcf().canvas.draw()
      except ValueError:
        pass
      finally:
        self.lines=[]

  def toDict(self):
    return {'Left': self.left, 'Right': self.right, 
        'Bottom': self.bottom, 'Top': self.top }

  def toFile(self, filename=None, mode='w'):
    self.filename = filename or self.filename
    FileIO.savesettings(self.filename, mode,
      **{self.name: self.toDict()} )

  def toAbsolute(self,origin):
    if self.origin != 'absolute':
      self.left+=origin[0]
      self.right+=origin[0]
      self.bottom+=origin[1]
      self.top+=origin[1]
      self.origin='absolute'
    return self

  def toRelative(self,origin):
    if self.origin == 'relative':
        return self
    elif self.origin == 'absolute':
      self.left-=origin[0]
      self.right-=origin[0]
      self.bottom-=origin[1]
      self.top-=origin[1]
      self.origin = 'relative'
      return self
    else:
        raise ROIError, "Origin must be either 'relative' or 'absolute'"
        
  @property
  def width(self):
    return self.right-self.left+1

  @property
  def height(self):
    return self.top-self.bottom+1

  def __repr__(self):
    if None in (self.left, self.right, self.bottom, self.top):
        name = self.name or 'Undefined'
        return "<ROI '%s' = uninitialized>" % self.name
    else:
        return "<%s ROI '%s' = L: %d, R: %d, B: %d, T: %d>" % \
      (self.origin,self.name,self.left,self.right,self.bottom,self.top) 


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

    self._showROI = False

    #################################################
    ## Load data from file called filename if string
    #################################################
    if isinstance(filename, str):
      self.filename = filename
      self._img = FileIO.loadimg(filename)
      self._roi = {}
      self._donorROIName = Stack.defaultDonorROI
      self._acceptorROIName = Stack.defaultAcceptorROI


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

      self.origin = (self.roileft,self.roibottom)
      self.addROI(*self.__class__.defaultROI.values())

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
      self.origin = filename.origin
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
  def time(self):
      return np.arange(1,self.frames+1)*self.exposurems/1000.

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
        roi = ROI.copy(roi)
        key = roi.name

        # recast to relative origin
        roi = roi.toRelative( self.origin )

        if roi.right > self.width or roi.top > self.height:
          raise ROIError, "ROI is larger than Stack dimensions"

        self._roi[key] = roi

      except AttributeError:
        raise TypeError, "Must use objects with ROI interface"
      except ROIError:
        raise ROIError, "ROI out of bounds! %s" % repr(roi)

  def showROI(self,*args):
    self._showROI=True
    for roi in args:
      self._roi[roi].draw()

  def show(self, frame, **kwargs):
    return self[frame].show(**kwargs)
 
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
      if self._roi.has_key(str(roi)):
        roi = self._roi[roi]
      return self[:,roi.bottom:roi.top,roi.left:roi.right].counts()
    else:
      return np.sum( np.sum(self._img,axis=1), axis=1 )

  def attime(self,time):
      if isinstance(time,slice):
        start,step = None,None
      if time.start:
        start = time.start/self.exposurems
      if time.step:
        step = time.step/self.exposurems
      time = slice(start,time.stop/self.exposurems,step)
      return self[time/self.exposurems]
      
  def __getitem__(self,key):
    if isinstance(key,int):
      return Frame(self._img[key], self._roi)
    else:
      temp = Stack(self)
      temp._img = temp._img[key]
      if isinstance(temp._img.shape,tuple) and len(temp._img.shape) > 2:
        temp.frames = temp._img.shape[0]
      else:
        temp.frames = 1

      return temp
    raise IndexError, "Invalid index: %s" % str(key)

  def append(self, stack):
    temp = copy.copy(self)
    temp._img = np.append( temp._img, stack._img, axis=0 )
    return temp

  def __sub__(self, stack):
    return self.__add__(stack.__neg__())

  def __add__(self, stack):
    temp = Stack(self)
    try:
      temp._img = temp._img + stack._img
    except ValueError:
      raise StackError("Couldn't add images: check sizes are the same")
    return temp

  def __neg__(self):
    temp = Stack(self)
    temp._img = -temp._img
    return temp

  def __repr__(self):
    return "Stack %dx%dx%d" % (self.frames, self.height, self.width)

  def __iter__(self):
    for i in range(self.frames):
      yield self[i]
    
class Frame:
    
  def __init__(self, imgarray, roi=None):
    self._img = imgarray
    self._roi = roi

  def __getitem__(self,key):
    try:
        return self._img[ key.bottom:key.top, key.left:key.right ]
    except AttributeError:
        return self._img[:,key]

  @property
  def height(self):
    return self._img.shape[0]

  @property
  def width(self):
    return self._img.shape[1]

  def counts(self,roi):
    return np.sum( self[roi] )

  def show(self, **kwargs):
    cmap = kwargs.get('cmap')

    plt.hold(1)
    plt.cla()
    p = plt.imshow(self._img, cmap=cmap)
    p.get_axes().invert_yaxis()
    if self._roi is not None:
      for roi in self._roi.itervalues():
        roi.draw()
    plt.draw()

  def __neg__(self):
    temp = copy.copy(self)
    temp._img = -temp._img
    return temp

  def __repr__(self):
    return "Image %dx%d" % self._img.shape

