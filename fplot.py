import os.path as path
from itertools import cycle
import matplotlib.pyplot as plt
from collections import OrderedDict, namedtuple

from datatypes import TrapData,hasTrapData,hasFretData
import constants

class Figure(object):
  def __init__(self, fig_id=None):
    self.figure_id = fig_id
    self._figure = None

  def __getstate__(self):
    ''' Return __dict__ for pickling with _figure attribute removed (can't be pickled)
    '''
    state = self.__dict__.copy()
    state['_figure'] = None
    return state

  @classmethod
  def fromCurrent(cls):
    gcf = plt.gcf()
    current = cls(gcf.number)
    current._figure = gcf
    return current

  def show(self):
    self.exists or self.new()
    plt.figure(self.figure_id)
    return self

  def new(self):
    self._figure = plt.figure(self.figure_id)
    self.figure_id = self._figure.get_label() or self._figure.number
    return self

  @property
  def visible(self):
    return self.exists and plt.fignum_exists(self._figure.number)

  @property
  def exists(self):
    return self._figure is not None

  @property
  def figure(self):
    if not self.exists:
      self.new()
    return self

  def pickPoints(self, num_of_pts=2):
  	return self._figure.ginput(num_of_pts)

  def pickRegions(self, num=1):
    points = sorted(x for x,f in self.pickPoints(num*2))
    return [(points[i],points[i+1]) for i in range(0,len(points),2)]
	
  def makeCurrent(self):
    if not self.exists:
      raise RuntimeError('Figure is not visible')
    plt.figure(self._figure.number)
    return self

  def plot(self, *args, **kwargs):
    self.show()
    try:
      # Treat the first argument as an object that can plot itself...
      return args[0].plot(*args[1:], **kwargs)
    except AttributeError:
      # ...unless it can't
      return plot(*args, **kwargs)

  def plotall(self, *args, **kwargs):   
    plotall(*args, **kwargs)

  def clear(self):
    if self.visible:
      self._figure.clf()
      plt.draw()
    return self
      
  def close(self):
    if self.exists:
      plt.close(self._figure)
    return self

  def annotate(self, text, location):
    "Annotate figure with text at location (x,y)"
    x,y = location
    return plt.text(x, y, text)

  IMAGE_OUTPUT_FORMATS = ('emf', 'eps', 'pdf', 'png', 'ps',
      'raw', 'rgba', 'svg', 'svgz') 

  DEFAULT_FILE_DIMENSIONS = (9, 7.5)

  def toFile(self, filename=None, size=None):
    size = size or Figure.DEFAULT_FILE_DIMENSIONS
    if filename:
      ext = path.splitext(filename)[1]
      if ext[1:] not in Figure.IMAGE_OUTPUT_FORMATS:
        filename += constants.DEFAULT_FIGURE_EXT
    else:
      filename = 'Figure {0}{1}'.format(self.figure.number, constants.DEFAULT_FIGURE_EXT)
    self._figure.set_size_inches(*size)
    self._figure.savefig(filename, bbox_inches='tight', pad_inches=0.1)

Region = namedtuple('Region', 'start end')

class Regions(list):
  def __init__(self, iterable):
      items = map(Region._make, iterable)
      if not self._valid_regions(items):
        raise ValueError('Regions must be non-overlapping: {}'.format(items))
      super(Regions, self).__init__(items)
      
  def _valid_regions(self, regions):
    last = 0
    for start,end in regions:
      if start < last:
        return False
    return True

  def __str__(self):
      return '\n'.join(
    'Region {n}: ({start:0.1f}, {end:0.1f})'.format(
      n=n, start=region.start, end=region.end) for n,region in enumerate(self)
    )

  def __add__(self, other):
      return Regions(Region(x+other, y+other) for x,y in self)
  
  def __setitem__(self, index, val):
      super(Regions, self).__setitem__(index, Region._make(val))


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

def _has_color(style_string):
  return style_string and len(style_string) > 0 and style_string[0].isalnum()

class PlotStyle(dict):
  """Dictionary-like class for setting styles on fplots

PlotStyle constructor will only allow setting styles on things that
fplot.plot() is going to plot: check PlotStyle.ITEMS

    style = PlotStyle(donor='r--', fret='.')
    fplot.plot(fret_data_object, style=style)

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
      color = '' # lets matplotlib decide
    default_style = [color+linestyle for linestyle in PlotStyle.STYLES]
    for item,style in styles.items():
      if style and not _has_color(style):
        styles[item] = color + style
    super(PlotStyle, self).__init__(
      zip(PlotStyle.ITEMS, default_style),
      **styles)

  @classmethod
  def with_default_style(cls, style):
    style_dict = dict(zip(cls.ITEMS, [style]*len(cls.ITEMS)))
    return cls(**style_dict)
    
  def __setitem__(self, key, value):
    if key not in PlotStyle.ITEMS:
      raise ValueError('Item "{}" is not plotted and has no style')
    super(PlotStyle, self).__setitem__(key, value)

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
  kwargs.setdefault('markersize', 2)

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
  donor, acceptor = 'donor', 'acceptor'
  if hasFretData(data) and displayFRET:
    plt.subplot(*next(layout))
    not hold and plt.cla()
    plt.hold(True)
    ax1 = subplot(data.time, data.donor, style[donor], label=donor, **kwargs)
    subplot(data.time, data.acceptor, style[acceptor], label=acceptor, axes=('','counts'), **kwargs)
    plt.hold(hold)
    if loc is not None:
      plt.legend(loc=loc,ncol=2,prop={'size':'small'})
  if displayFRET:
    subplot(data.time, data.fret, style['fret'], layout=next(layout), 
              axes=('Seconds','FRET'))
    plt.ylim(-0.1, 1.1)

  ax2 = None
  if pull:
    trap_style = PlotStyle.TRAP_STYLE_NAME
    x_coord,x_label = (pull.ext,'Extension (nm)') if FEC else (pull.sep,'Separation (nm)')
    ax2 = subplot(x_coord, pull.f, style[trap_style], layout=next(layout), axes=(x_label,'Force (pN)'), label=label, **kwargs)
    if loc is not None:
      plt.legend(loc=loc,ncol=2,prop={'size':'small'})

  first_plot = ax1 or ax2
  first_plot.set_title(title)

fret_default_style = OrderedDict(
                    [('donor', '--'), 
                    ('acceptor', '-'),
                    ('fret', '-')])

def _num_of_visible_styles(style_dict):
  return len(filter(None, style_dict.values()))

def _n_row_layout_generator(num):
  return ((num,1,x) for x in range(1,num+1))

def axes_label(x, y):
  plt.xlabel(x)
  plt.ylabel(y)

def axes_limits(x=None, y=None):
  return plt.xlim(x), plt.ylim(y)
  
def plot_fret(data, hold=False, layout=None, **styles):
  styles_w_defaults = fret_default_style.copy()
  styles_w_defaults.update(styles)
  donor, acceptor, fret = styles_w_defaults.values()
  if layout is None:
    num_layouts = 2 if fret and (donor or acceptor) else 1
    layout = _n_row_layout_generator(num_layouts)
  elif isinstance(layout, list):
    layout = iter(layout)
  else:
    layout = cycle([layout])

  not hold and plt.cla()
  plt.hold(True)

  count_layout = next(layout)
  if donor:
    subplot(data.time, data.donor, donor, label='donor', layout=count_layout)
    count_layout = None

  if acceptor:
    subplot(data.time, data.acceptor, acceptor, label='acceptor', layout=count_layout)
  axes_label('', 'Counts')

  if fret:
    subplot(data.time, data.fret, fret, layout=next(layout), 
            axes=('Seconds','FRET'))
    plt.ylim(-0.1, 1.1)
   
  plt.legend(loc='best', ncol=2, prop={'size':'small'})

def subplotsNeeded(data):
  num = 0
  if hasattr(data, 'donor') and hasattr(data, 'acceptor'):
    num += 1
  if hasattr(data,'fret'):
    num += 1
  if hasTrapData(data):
    num += 1
  return num

def subplot(*args, **kwargs):
  layout = kwargs.pop('layout',())
  axes = kwargs.pop('axes', ())
  hold = kwargs.get('hold', None)

  if hold is not None:
    plt.hold(hold)
  if layout:
    ax = plt.subplot(*layout)
  else:
    ax = plt.gca()
  plt.plot(*args,**kwargs)
  ax.autoscale_view(tight=True)
  if axes:
    try:
      plt.xlabel(axes[0])
      plt.ylabel(axes[1])
    except IndexError:
      raise ValueError('subplot expects labels for BOTH axes')
  return ax

def hist(data, bins=50, plot='fret', hold=False, **kwargs):
  plt.hold(hold)
  counts,bins,patches = plt.hist( getattr(data, plot), bins )
  delta = bins[1]-bins[0]
  bins = (bins-delta/2)[1:]
  return bins,counts
