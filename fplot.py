import os.path as path
from itertools import cycle
from collections import defaultdict
import matplotlib.pyplot as plt

from datatypes import TrapData,hasTrapData,hasFretData,hasTrapData
import constants

class Figure(object):
  def __init__(self, fignumber=None):
    self.figure = plt.figure(fignumber) if fignumber else None

  @classmethod
  def fromCurrent(cls):
    return cls(plt.gcf().number)

  def new(self):
    self.figure = plt.figure()

  @property
  def exists(self):
    return self.figure is not None and plt.fignum_exists(self.figure.number)

  def pickPoints(self, num_of_pts=2):
  	return self.figure.ginput(num_of_pts)

  def pickRegions(self, num=1):
    points = sorted(x for x,f in self.pickPoints(num))
    return [(points[i],points[i+1]) for i in range(0,len(points),2)]
	
  def makeCurrent(self):
    if not self.exists:
      raise RuntimeError('Figure object does not exist')
    plt.figure(self.figure.number)
    return self

  def plot(self, *args, **kwargs):
    if not self.exists:
      self.new()
    else:
      self.makeCurrent()
    try:
      # Treat the first argument as an object that can plot itself...
      return args[0].plot(*args[1:], **kwargs)
    except AttributeError:
      # ...unless it can't
      return plot(*args, **kwargs)

  def clear(self):
    if self.exists:
      self.figure.clf()
      self.figure.show()
      
  def close(self):
    if self.exists:
      plt.close(self.figure)

  def annotate(self, text, location):
    "Annotate figure with text at location (x,y)"
    x,y = location
    return plt.text(x, y, text)
  IMAGE_OUTPUT_FORMATS = ('emf', 'eps', 'pdf', 'png', 'ps',
      'raw', 'rgba', 'svg', 'svgz') 

  DEFAULT_SIZE = (9, 7.5)
  def toFile(self, filename=None, size=None):
    size = size or Figure.DEFAULT_SIZE
    if filename:
      ext = path.splitext(filename)[1]
      if ext[1:] not in Figure.IMAGE_OUTPUT_FORMATS:
        filename += constants.DEFAULT_FIGURE_EXT
    else:
      filename = 'Figure {0}{1}'.format(self.figure.number, constants.DEFAULT_FIGURE_EXT)
    self.figure.set_size_inches(*size)
    self.figure.savefig(filename, bbox_inches='tight', pad_inches=0.1)

def plotall(fret, pull=None,  **kwargs):
  assert isinstance(fret, (list, tuple, type(None)))
  if pull is not None and len(fret) != len(pull):
    raise ValueError('Number of pull ({}) must match number of FRET objects')
  elif pull is None:
    pull = [None]*len(fret)

  labels = kwargs.pop('labels', [])
  for obj in fret:
    label = labels.pop(0) if len(labels)>0 else ''
    plot(obj, pull=pull.pop(0), label=label, **kwargs)


class PlotStyle(object):
  """Dictionary-like class for setting styles on fplots

PlotStyle constructor will only allow setting styles on things that
fplot.plot() is going to plot: check PlotStyle.ITEMS

>>> style = PlotStyle(donor='r--', fret='.')
>>> fplot.plot(fret_data_object, style=style)

    # donor is dots with default color scheme
    plot( ... , ... , style=PlotStyle(donor='.') )
    
    # donor is blue line, acceptor is red line, rest is black
    plot( ... , ... , style=PlotStyle(color='k', donor='b-', acceptor='r-') )
    
    # all plots are lines
    plot( ... , ... , style='-' )
    
    # all plots are green
    plot( ... , ... , style=PlotStyle(color='g') )
    
    # color cycles for every line drawn (like matplotlib default)
    plot( ... , ... , style=PlotStyle(color='auto') )
"""

  ITEMS = ('donor', 'acceptor', 'fret', 'trap')
  STYLES = (':','-','-','.')

  TRAP_STYLE_NAME = 'trap'

  def __init__(self, color=None, **styles):
    if color is None:
      color = next_color_in_cycle()
    elif color == 'auto':
      color = ''
    default_style = [color+linestyle for linestyle in PlotStyle.STYLES]
    for item,style in styles.items():
      if not style[0].isalnum():
        styles[item] = color+style
    self._styles = dict(
      zip(PlotStyle.ITEMS, default_style),
      **styles)

  @classmethod
  def with_default_style(cls, style):
    style_dict = dict(zip(cls.ITEMS, [style]*len(cls.ITEMS)))
    return cls(**style_dict)

  def __getitem__(self, key):
    return self._styles[key]

  def __setitem__(self, key, value):
    if key not in PlotStyle.ITEMS:
      raise ValueError('Item "{}" is not plotted and has no style')
    self._styles[key] = value


COLOR_CYCLE = plt.rcParams['axes.color_cycle']
COLOR = (color for color in cycle(COLOR_CYCLE))
next_color_in_cycle = lambda : next(COLOR)

def plot(data, pull=None, style=None, **kwargs):
  """Plot FretData and/or TrapData as stacked subplots of counts, FRET, and FEC/FDC
  @data: datatypes.AbstractData
  @pull: datatypes.AbstractData
  @style: dict or str
  """
  if isinstance(style, str):
    style = PlotStyle.with_default_style(style)
  elif style is None:
    style = PlotStyle()
  loc = kwargs.pop('legend', 'best')
  title = kwargs.pop('title','')
  label = kwargs.pop('label', '')
  displayFRET = kwargs.pop('show_fret', hasattr(data,'fret'))

  hold=kwargs.pop('hold', None)
  if hold is not None:
    plt.hold(hold)

  if not pull and hasTrapData(data):
    pull = TrapData.fromObject(data)

  num = 0
  if displayFRET and hasattr(data, 'donor'): num += 1
  if displayFRET and hasFretData(data): num += 1
  if pull:
    num += 1
    FEC = kwargs.pop('FEC', num==1)
  else:
    FEC = kwargs.pop('FEC', False)

  if num == 0:
    raise ValueError("Don't know how to plot arguments: need TrapData or FretData")

  if data is None and pull is not None:
    layout = iter([(3,1,3)])
  else:
    layout = iter((num,1,x) for x in range(1,num+1))

  ax1 = None
  if hasFretData(data) and displayFRET:
    donor, acceptor = 'donor', 'acceptor'
    plt.subplot(*next(layout))
    not hold and plt.cla()
    plt.hold(True)
    ax1 = _subplot(data.time, data.donor, style[donor], label=donor, **kwargs)
    _subplot(data.time, data.acceptor, style[acceptor], label=acceptor, axes=('','counts'), **kwargs)
    plt.hold(hold)
    if loc is not None:
      plt.legend(loc=loc,ncol=2,prop={'size':'small'})
    if displayFRET:
      _subplot(data.time, data.fret, style['fret'], layout=next(layout), 
                axes=('Seconds','FRET'))
      plt.ylim(-0.1, 1.1)

  ax2 = None
  if pull:
    trap_style = PlotStyle.TRAP_STYLE_NAME
    x_coord,x_label = (pull.ext,'Extension (nm)') if FEC else (pull.sep,'Separation (nm)')
    ax2 = _subplot(x_coord, pull.f, style[trap_style], layout=next(layout), axes=(x_label,'Force (pN)'), label=label, **kwargs)
    if loc is not None:
      plt.legend(loc=loc,ncol=2,prop={'size':'small'})

  first_plot = ax1 or ax2
  first_plot.set_title(title)

def subplotsNeeded(data):
  num = 0
  if hasattr(data, 'donor') and hasattr(data, 'acceptor'):
    num += 1
  if hasattr(data,'fret'):
    num += 1
  if hasTrapData(data):
    num += 1
  return num

def _subplot(*args, **kwargs):
  sub = kwargs.pop('layout',())
  axes = kwargs.pop('axes', ())
  hold = kwargs.get('hold', None)

  if hold is not None:
    plt.hold(hold)
  if sub:
    ax = plt.subplot(*sub)
  else:
    ax = plt.gca()
  plt.plot(*args,**kwargs)
  ax.autoscale_view(tight=True)
  if axes:
    try:
      plt.xlabel(axes[0])
      plt.ylabel(axes[1])
    except IndexError:
      raise ValueError('_subplot expects labels for BOTH axes')
  return ax

def hist(data, bins=50, plot='fret', hold=False, **kwargs):
  plt.hold(hold)
  counts,bins,patches = plt.hist( getattr(data, plot), bins )
  delta = bins[1]-bins[0]
  bins = (bins-delta/2)[1:]
  return bins,counts
