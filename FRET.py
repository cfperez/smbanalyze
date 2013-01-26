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
    x_coord = data.ext if FEC else data.sep
    _subplot(x_coord, data.f, '.', layout=next(layout), axes=('Sep (nm)','Force (pN)'))

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

def oldplot(*data, **kwargs):
  """loc=legend location (or None for off), title, hold (False for individual plots),
  names=labels for individual traces (must equal number of traces)
  prefix=for all traces in legend"""

  loc = kwargs.get('legend', 'best')
  hold = kwargs.get('hold',True)

  if len(data) == 1 and isinstance(data[0],list):
    data = data[0]

  # to prefix labels in plots
  names = kwargs.get('names',[None]*len(data))
  prefix = kwargs.get('prefix','')
  title = kwargs.get('title','')

  if len(names) != len(data):
    raise ValueError, "Must have same number of names as traces to plot"

  for name,trace in zip(names,data):
    hasfret = kwargs.get('fret',hasattr(trace,'fret'))
    if kwargs.get('FDC', hasattr(trace,'f')):
      FDC = (trace.sep, trace.f)
    else:
      FDC = None

    if not hold:
      plt.figure()

    fig_size = 1
    if hasfret is not False:
      fig_size += 1
    if FDC:
      fig_size += 1

    #if hasfret is not False:
    current_panel = 1
    s = plt.subplot(fig_size,1,current_panel)
    s.autoscale_view(tight=True)
    plt.title(title)

    name = name or getattr(trace,'molname','')

    if FDC:
      plt.plot(*FDC,label=prefix+' '+name)
      current_panel += 1
      s=plt.subplot(fig_size,1,current_panel)
      s.autoscale_view(tight=True)

    x_axis = None
    if hasattr(trace,'time'):
      x_axis = trace.time
      plt.xlabel('Seconds')
    elif not hasattr(trace, 'donor'):
      continue
    else:
      x_axis = range(1,len(trace.donor)+1)
      plt.xlabel('Frames')

    plt.plot(x_axis, trace.donor, label=' '.join([prefix,name,'donor']))
    plt.plot(x_axis, trace.acceptor,label=' '.join([prefix,name,'acceptor']))
    plt.ylabel('Counts')
    if loc:
      plt.legend(loc=loc,ncol=2,prop={'size':'small'})

    if hasfret is not False:
      current_panel += 1
      s=plt.subplot(fig_size,1,current_panel)
      s.autoscale_view(tight=True)
      plt.ylabel('FRET')
      plt.xlabel('Frames')
      if hasfret is True:
        plt.plot( trace.fret, label=' '.join([prefix,name]))
      else:
        plt.plot(hasfret, label=' '.join([prefix,name]))

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

def calculate(stack, beta=Constants.beta, gamma=Constants.gamma, minsub=False):
  """Calculates FRET of a pull from an Image.Stack

  calculate( Image.Stack, beta = Image.beta, gamma = Image.gamma)

  RETURNS array of calculated FRET for each frame
  """

  donor = stack.donor - (minsub and min(stack.donor))
  acceptor = stack.acceptor - donor*beta
  acceptor = acceptor - (minsub and min(acceptor))

  return FretData(stack.time, donor, acceptor, acceptor/(acceptor+gamma*donor))

def calcToFile(stack, filename, **kwargs):
  "saveFRETdata(fret, ImageStack, filename): saves donor,acceptor, FRET to 3 column text file"

  fretdata = calculate(stack, **kwargs)
  FileIO.savedat(filename, (stack.time,stack.donor,stack.acceptor,fretdata), header='time donor acceptor FRET', fmt=('%.3f','%u','%u','%.5f'))
  return fretdata

def toFile(filename, data):
  try:
    FileIO.savedat(filename, (data.time,data.donor,data.acceptor,data.fret), header='time donor acceptor FRET', fmt=('%.3f','%u','%u','%.5f'))
  except AttributeError:
    raise AttributeError('FRET.save expects argument with time, donor, acceptor, and fret attributes')

def fromFile(filename, **kwargs):
  return FileIO.loadfret(filename, **kwargs)

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
    toFile(FileIO.change_extension(fname,ext), output)
    all_output += [output]

  return all_output

def calcDirectory(select='',*args, **kwargs):

  dir = kwargs.get('dir')
  roi_origin = kwargs.get('roi_origin','absolute')
  verbose = 'verbose' in args or kwargs.get('verbose')
  plotall = 'plotall' in args or kwargs.get('plotall')
  hold = 'hold' in args or kwargs.get('hold')
  saveplot = 'saveplot' in args or kwargs.get('saveplot')
  aslist = 'aslist' in args or kwargs.get('aslist',False)
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
        *Image.ROI.fromFile(roi_file, origin=roi_origin))
      if verbose:
        print "Using ROI file: %s" % roi_file

    if not Image.Stack.defaultROI:
      raise RuntimeError, "No ROIs set or loaded--cannot compute counts"

    img_files = glob.glob('*.img')

    if verbose:
      print "Found files:\n%s\n" % '\n'.join(img_files)

    BG = None             # current background image
    results = [] if aslist else useful.dotdict() # output

    for fname,background in matchImgFilestoBackground(img_files):

      if select not in fname:
        continue

      basename, ext = os.path.splitext(fname)

# Get molecule information from filename
      finfo = FileIO.parseFilename(basename)
      (construct, context, slide, mol, pull, force, min, sec,
  series, isBackground) = finfo

# Only process requested slide and molecules
      if (not match_or_included(user_slide,slide) or 
        not match_or_included(user_mol,mol)):
        if verbose:
          print "skipping " + fname
        continue


      image = Image.Stack(fname)

      if background:
        # Load the background file from disk only if needed
        if BG is None or background != BG.filename:
          BG = Image.fromFile(background, background=True)
        image = image - BG

      if verbose:
        print "Processing...\n\timage: %s\n\tbackground: %s\n" \
              % (fname,background)

      data = calcToFile( image, basename+'.fret' )

      temp = useful.dotdict(image=image, fret=data, 
        donor=image.donor, acceptor=image.acceptor,molname=molname(finfo))

      if aslist:
        temp.info = finfo
        results += [temp]
      else:
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

  if not aslist:
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


def matchImgFilestoBackground(img_files=None):
  "Look through current directory and intelligently match each img file to a background.img file"

  bg_files = glob.glob('*_background.img')
  img_files = img_files or glob.glob('*.img')

  background = ''
  slide_background = '' # last background used on this slide
  last_slide = None     # slide number of last file

  matched_files = []
  bg = []

  for fname in img_files:

	if '_background' in fname:
	  continue

	basename, ext = os.path.splitext(fname)

	# Get molecule information from filename
	finfo = FileIO.parseFilename(fname)
	(construct, context, slide, mol, pull, force, min, sec,
		series, isBackground) = finfo

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

	last_slide = slide

	bg += [background]
	matched_files += [fname]

  return zip(matched_files,bg)
