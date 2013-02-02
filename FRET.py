import os
import glob
import re
import operator
from itertools import izip

import matplotlib.pyplot as plt
from numpy import concatenate

import useful
import Image
import FileIO
import Constants
from Types import PullFretData,FretData,hasPullData,hasFretData

molID = lambda t: 's{0}m{1}'.format(t.slide,t.mol)
molname = lambda t: 's{0}m{1}_{2}'.format(t.slide,t.mol,t.pull)
pN = lambda f: 'f'+str(f)+'pN'

def processFiles(flist, roi='roi.txt', background=None, ext=FileIO.FRET_FILE):
  "processFiles(filelist, roi='roi.txt', background=None, ext=Fret_File_extension)"

  if background:
    BG = Image.fromFile(background,background=True)
  else:
    BG = Constants.default_background_subtract

  if isinstance(roi,str):
    roi = Image.ROI.fromFile(roi)

  all_output = []
  for fname in flist:
    img = Image.fromFile(fname) - BG
    img.addROI(*roi)
    output = calculate(img)
    toFile(FileIO.change_extension(fname,ext), output, img.metadata)
    all_output += [output]

  return all_output

def multiplot(*data, **kwargs):
  names = kwargs.get('names',[None]*len(data))
  prefix = kwargs.get('prefix','')
  title = kwargs.get('title','')

  if not isinstance(datalist,list):
    datalist=[datalist]

  if not (isinstance(titles,tuple) or isinstance(titles,list)):
    titles = (titles,)

  if len(titles) < len(datalist):
    titles += (titles[-1],) * (len(datalist)-len(titles))

  figures=[]
  for data,title in izip(datalist,titles):
    if len(datalist)>1 or kwargs.get('autofigure', False):
      figures += [plt.figure()]
    elif not kwargs.get('hold',False):
      plt.clf()
      figures += [plt.gcf()]

  return figures

def plot(data, pull=None, **kwargs):
  loc = kwargs.get('legend', 'best')
  title = kwargs.get('title','')
  FEC = kwargs.get('FEC',False)

  hold=kwargs.get('hold',False)
  plt.hold(hold)

  if pull and not hasPullData(data):
    data = PullFretData(*(pull+data))

  num = kwargs.get('numplot',subplotsNeeded(data))
  if num==0:
    raise ValueError("Don't know how to plot argument: missing named fields")

  layout = iter((num,1,x) for x in range(1,num+1))

  if hasFretData(data):
    plt.subplot(*next(layout))
    not hold and plt.cla()
    plt.hold(True)
    _subplot(data.time, data.donor, label='donor')
    _subplot(data.time, data.acceptor, label='acceptor',axes=('','counts'))
    plt.hold(hold)
    plt.title(title)
    plt.legend(loc=loc,ncol=2,prop={'size':'small'})

  if hasattr(data,'fret'):
    _subplot(data.time, data.fret, layout=next(layout), axes=('Seconds','FRET'))

  if hasPullData(data):
    x_coord,x_label = (data.ext,'Extension (nm)') if FEC else (data.sep,'Separation (nm)')
    _subplot(x_coord, data.f, '.', layout=next(layout), axes=(x_label,'Force (pN)'))

def subplotsNeeded(data):
  num = 0
  if hasFretData(data):
    num += 2
  elif hasattr(data,'fret'):
    num += 1
  if hasPullData(data):
    num += 1

  return num

def _subplot(*args,**kwargs):
  sub = kwargs.pop('layout',())
  axes = kwargs.pop('axes', ())

  if sub:
    plt.subplot(*sub)
  plt.plot(*args,**kwargs)
  plt.gca().autoscale_view(tight=True)
  if axes:
    try:
      plt.xlabel(axes[0])
      plt.ylabel(axes[1])
    except IndexError:
      raise ValueError('_subplot expects labels for BOTH axes')

def hist(*data, **kwargs):
  bins = kwargs.get('bins',50)
  attr = kwargs.get('plot','fret')

  if kwargs.get('join',True):
    func = concatenate
  else:
    func = lambda x: x

  counts,bins,patches = plt.hist(
        func(map(operator.attrgetter(attr),data)), bins)

  delta = bins[1]-bins[0]
  bins = (bins-delta/2)[1:]

  return bins,counts

def calcToFile(stack, filename, **kwargs):
  "Calculates and saves saves donor, acceptor, calculated FRET values to 3 column text file"
  fretdata = calculate(stack, **kwargs)
  toFile(filename, fretdata, stack.metadata)
  return fretdata

def calculate(stack, beta=Constants.beta, gamma=Constants.gamma, minsub=False):
  """Calculates FRET of a pull from an Image.Stack

  calculate( Image.Stack, beta = Constants.beta, gamma = Constants.gamma)

  RETURNS array of calculated FRET for each frame
  """
  donor = stack.donor - (minsub and min(stack.donor))
  acceptor = stack.acceptor - donor*beta
  acceptor = acceptor - (minsub and min(acceptor))
  return FretData(stack.time, donor, acceptor, acceptor/(acceptor+gamma*donor))

def toFile(filename, data, metadata, comments=''):
  return FileIO.savefret(filename, data, metadata, comments)

def fromFile(filename, **kwargs):
  return FileIO.load(filename, **kwargs) #comments=FileIO.toSettings, **kwargs)
