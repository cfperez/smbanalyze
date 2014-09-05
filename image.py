from __future__ import with_statement
from operator import methodcaller,itemgetter,add as opadd
from collections import defaultdict
from itertools import cycle, dropwhile
import copy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import fileIO
from figure import Figure
import useful
import constants

################################################################################
## EXCEPTIONS
################################################################################
class StackError(Exception):
  pass

class ROIError(Exception):
  pass

class ROIOutofBounds(ROIError):
  pass

class ROIInvalid(ROIError):
  pass

def setDefaultROI(*args):
  Stack.setDefaultROI(*args)

def fromFile(filename, background=None, **kwargs):
  """fromFile(filename, [background='' or True]): 
  
  Load Image from filename with optional background file or constant to subtract.

  fromFile('image.img', background = 'filename') => load background file and subtract (SLOW!)
  fromFile('image.img', background=100) => subtract 100 from image as background level
""" 
  roi = kwargs.pop('roi', None)
  img = Stack.fromfile(fileIO.add_img_ext(filename))

  bgimg = fromBackground(background) if background else ConstantStack(0)
  bg_exposure = bgimg.metadata['exposurems']
  if not isinstance(bgimg, ConstantStack) and bg_exposure != img.metadata['exposurems']:
    print '''WARNING Trying to subtract background with exposure time {} 
      from image with exposure time {}'''.format(
        bg_exposure,
        img.metadata['exposurems'])
  img.metadata['background'] = bgimg.metadata['filename'] or background
  img -= bgimg

  if roi is not None:
    if isinstance(roi, basestring):
      roi = ROI.fromFile(roi)
    img.addROI(*roi)
  return img

def is_background(bg):
  pass

def fromBackground(name=None, zfilter='median'):
  if isinstance(name, str):
    bg = Stack.fromfile(name)
  elif isinstance(name, Stack):
    bg = name
  elif isinstance(name, int) or name is None:
    return ConstantStack(name or constants.default_background_subtract)
  else:
    raise ValueError('Argument "name" must be either a string, stack, int, or None')
  return bg.toBackground(zfilter)
    

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

  def __init__(self, bottomLeft, topRight, name='', origin='relative'):
    L, B = map(int, bottomLeft)
    R, T = map(int, topRight)
    if R <= L:
      raise ROIInvalid('Right {0} must be greater than left {1}'.format(R,L))
    if T <= B:
      raise ROIInvalid('Top {0} must be greater than bottom {1}'.format(T,B))
    if R<0 or L<0 or T<0 or B<0:
      msg = 'Boundaries must always be >0: {0}'.format((L,R,T,B))
      raise ROIOutofBounds(msg)
    self.left = L
    self.right = R
    self.top = T
    self.bottom = B
    self.name = name
    self.origin = origin
    self.lines = []

  @classmethod
  def fromFile(cls, filename, origin=None):
    "Load ROI(s) from a LabView config file and return tuple of Image.ROI objects"
    settings = fileIO.loadsettings(filename)
    temp = []
    origin = origin or settings.pop('origin', 'absolute')
    for name, roi in settings.items():
      corners = itemgetter('left','right','bottom','top')
      corners_int = map(useful.toInt, corners(roi))
      temp.append(
        cls.fromCorners(*corners_int, name=name, origin=origin)
      )
    return tuple(temp)

  @classmethod
  def copy(cls, roi, name='', origin=None):
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
  def fromDict(cls, d):
    return cls.fromCorners(d['left'],d['right'],d['bottom'],d['top'],
      name=d.get('name',''),origin=d.get('origin','relative'))

  @classmethod
  def save(cls, filename, *ROIs):
    "save(filename, *ROIs): Save ROI(s) to a LabView config file format"
    mode = 'w'
    for roi in ROIs:
      roi.toFile(filename, mode)
      mode = 'a'

  Colors = ('r','r')
  _lastColor = 0

  def draw(self, color=None):
    if not color:
      i = ROI._lastColor
      color=ROI.Colors[i]
      ROI._lastColor = (i+1) % (len(ROI.Colors))

    self.clear()

    self.lines = [
      plt.gca().add_patch( 
        Rectangle( (self.left,self.bottom), self.width, self.height, 
          fill=False, lw=1, color=color)
      )
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
    return {'left': self.left, 'right': self.right, 
        'bottom': self.bottom, 'top': self.top,
        'name': self.name, 'origin': self.origin }

  def toFile(self, filename=None, mode='w'):
    self.filename = filename or self.filename
    fileIO.savesettings(self.filename, mode,
      **{self.name: self.toDict()} )

  def _convert(self, origin):
    corners = map(opadd,
      (self.left,self.bottom,self.right,self.top), origin*2)
    return corners[0], corners[2], corners[1], corners[3]
    
  def toAbsolute(self, origin_coord):
    if self.origin != 'absolute':
      corners = self._convert(origin_coord)
      return ROI.fromCorners(*corners, name=self.name, origin='absolute')
    else: 
      return ROI.copy(self)

  def toRelative(self, origin_coord):
    if self.origin == 'relative':
        return self
    elif self.origin == 'absolute':
      try:
        corners = self._convert(useful.negate(origin_coord))
        return ROI.fromCorners(*corners, name=self.name, origin='relative')
      except ROIOutofBounds as e:
        if corners[0]<origin_coord[0]:
          raise ROIOutofBounds("ROI left ({0}) must be greater than subimage origin {1}".format(self.left, origin_coord[0]))
        elif corners[2]<origin_coord[1]:
          raise ROIOutofBounds("ROI right ({0}) must be greater than subimage origin {1}".format(self.bottom, origin_coord[1]))
        else:
          raise ROIOutofBounds("Check that ROI left and bottom is > origin:\n"+str(e))
    else:
        raise ROIError, "Origin must be either 'relative' or 'absolute'"
        
  @property
  def width(self):
    return self.right-self.left+1

  @property
  def height(self):
    return self.top-self.bottom+1

  @property
  def size(self):
    return self.height*self.width

  @property
  def edges(self):
    return self.left, self.right, self.bottom, self.top

  @property
  def corners(self):
    return ((x,y) for x in (self.left,self.right) for y in (self.bottom,self.top))

  def __repr__(self):
    if None in self.edges:
        name = self.name or 'Undefined'
        return "<ROI '%s' = uninitialized>" % name
    else:
        return ROI.REPR(self)

  REPR = lambda self: "<%s ROI '%s' = L: %d, R: %d, B: %d, T: %d>" % \
    (self.origin,self.name,self.left,self.right,self.bottom,self.top) 


class Stack:
  "A class to load and manipulate images generated during single molecule experiments"

  defaultROI = {}
  defaultDonorROI = 'donor'
  defaultAcceptorROI = 'acceptor'

  @classmethod
  def fromfile(cls, imgfile, camfile=''):
    img = fileIO.loadimg(imgfile)
    camfile = camfile or fileIO.to_cam_filename(imgfile)
    metadata = fileIO.loadcam(camfile)
    metadata['filename'] = imgfile
    stack = cls(img, metadata)
    # Keeping the .filename attribute for backward compatibility, but should move to
    # metadata lookup
    stack.filename = metadata['filename']
    return stack

  # def __init__(self, filename, camFile='', deepcopy=False):
  def __init__(self, data, metadata={}, roi={}):
    self._img = copy.copy(data)
    self.metadata = dict(metadata)
    self._donorROIName = Stack.defaultDonorROI
    self._acceptorROIName = Stack.defaultAcceptorROI
    self._figure = Figure()
    self._roi = roi
    if not roi:
      self.addROI(*self.defaultROI.values())
    self._frame_iter = cycle(range(self.frames))
    if metadata:
      self.origin = (self.metadata['roileft'],self.metadata['roibottom'])
      if self._img.shape != (self.metadata['frames'],self.metadata['height'],self.metadata['width']):
        raise StackError, ".img file and .cam file dimensions do not agree"

  def copy(self, deep=True):
    newstack = copy.copy(self)
    if deep:
      newstack._img = copy.copy(self._img)
    return newstack

  @property
  def frames(self):
    return self._img.shape[0]

  def __len__(self):
    return self.frames

  @property
  def time(self):
      return np.arange(1,self.frames+1)*self.metadata['exposurems']/1000.

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
    
  @property
  def roi(self):
    return self._roi

  @classmethod
  def setDefaultROI(cls, *args):
    for _roi in args:
        cls.defaultROI[_roi.name]=_roi

  def toBackground(self,zfilter='median'):
    width, height = self.width, self.height
    if zfilter == 'median':
      self._img = np.median( self._img, axis=0, overwrite_input=True ).reshape((1,height,width))
    elif zfilter == 'mean':
      self._img = np.mean( self._img, axis=0 ).reshape((1,height,width))
    elif zfilter == 'min':
      self._img = np.min( self._img, axis=0 ).reshape((1,height,width))
    else:
      raise ValueError, "Filter type can only be median, mean, or min"
    return self

  def addROI(self, *ROIs):
    for roi in ROIs:
      try:
        roi = ROI.copy(roi)
        key = roi.name
        roi = roi.toRelative(self.origin)
        if roi.right > self.width:
          raise StackError(
            "ROI 'right' {0} is outside right edge of image {1}: \n {2}".format(roi.right,self.width,roi)
          )
        if roi.top > self.height:
          raise StackError, "ROI 'top' is outside top edge of image: {0}\n {1}".format(roi.top,roi)
        self._roi[key] = roi
      except AttributeError:
        raise TypeError, "Must use objects with ROI interface"

  def showROI(self,*args):
    for roi in args:
      self._roi[roi].draw()

  def show(self, frame=None, **kwargs):
    if not isinstance(frame, int) and frame is not None:
      raise ValueError('First argument frame must be an integer')
    self._figure = kwargs.pop('figure', self._figure)
    if frame is None:
      frame = next(self._frame_iter)
    else:
      self._frame_iter = dropwhile(lambda n: n<=frame, cycle(range(self.frames)))
    self._figure.show()
    self._figure.makeCurrent()
    plt.title('Frame %d' % frame)
    self[frame].show(**kwargs)
    return frame
 
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
      roi = roi.toRelative(self.origin)
      return self[:,roi.bottom:roi.top,roi.left:roi.right].counts()
    else:
      return np.sum( np.sum(self._img,axis=1), axis=1 )

  def attime(self,time):
      if isinstance(time,slice):
        start,step = None,None
      exposurems = self.metadata['exposurems'] 
      if time.start:
        start = time.start/exposurems
      if time.step:
        step = time.step/exposurems
      time = slice(start,time.stop/exposurems,step)
      return self[time/exposurems]
      
  def __getitem__(self,key):
    if isinstance(key,int): # Single frame
      return Frame(self._img[key], self._roi)
    else: # It's a slice
      temp = self.copy(deep=False)
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
    temp = self.copy()
    if hasattr(stack,'_img'):
      try:
        temp._img = temp._img + stack._img
      except ValueError:
        raise StackError("Couldn't add images: check sizes are the same")
    else:
      temp._img = temp._img + stack
    return temp

  def __neg__(self):
    temp = self.copy()
    temp._img = -temp._img
    return temp

  def __eq__(self, other):
    return np.all(self._img == other._img)

  def __ne__(self, other):
    return not self==other

  def __repr__(self):
    return "Stack %dx%dx%d" % (self.frames, self.height, self.width)

  def __iter__(self):
    for i in range(self.frames):
      yield self[i]
    
class ConstantStack(Stack):
  def __init__(self, constant, metadata={}):
    Stack.__init__(self, constant, metadata)
    if not self.metadata:
      self.metadata = defaultdict(lambda : None)
    self._constant = constant

  @property
  def frames(self):
    return 1

  @property
  def time(self):
      return np.array([0])

  @property
  def height(self):
    return 0

  @property
  def width(self):
    return 0

  def __repr__(self):
    return "ConstantStack(%d)" % self._constant

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

  def counts(self, roi):
    return np.sum( self[roi] )

  def show(self, **kwargs):
    cmap = kwargs.get('cmap')
    plt.hold(1)
    plt.cla()
    p = plt.imshow(self._img, cmap=cmap)
    p.get_axes().invert_yaxis()
    if self._roi is not None and kwargs.get('roi', True):
      for roi in self._roi.itervalues():
        roi.draw()
    plt.draw()

  def __neg__(self):
    temp = copy.copy(self)
    temp._img = -temp._img
    return temp

  def __repr__(self):
    return "Image %dx%d" % self._img.shape