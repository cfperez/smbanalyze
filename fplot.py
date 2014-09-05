

__all__ = ['fec','fecs','fretfec','segmented','stack','stacks','fret','twinxaxis']

from itertools import cycle
import matplotlib.pyplot as plt
from collections import OrderedDict
from experiment import to_fret_ext
from datatypes import hasTrapData, TrapData

## Define some pretty colors
baby_blue = '#8ac7e6'
blue = '#17a0e6'
lightpink = '#f27474'
pink = '#f29d9d'
red = '#d92b2b'
salmon = '#f77254'
lightsalmon = '#f29f8d'
lightersalmon = '#fabdaf'
steel_blue = '#3296d2'
green = '#0b9e28'
lightgreen = '#2cde50'

COLOR_CYCLE = plt.rcParams['axes.color_cycle']
COLOR = (color for color in cycle(COLOR_CYCLE))
next_color_in_cycle = lambda : next(COLOR)

class Style(dict):
  """Dictionary-like class for setting styles on fplots

  Style constructor will only allow setting styles on things that
  fplot.plot() is going to plot: check Style.ITEMS

      style = Style(donor='r--', fret='.')
      fplot.plot(fret_data_object, style=style)

      # donor is dots with default color scheme
      plot( ... , ... , style=Style(donor='.') )
      
      # donor is blue line, acceptor is red line, rest is black
      plot( ... , ... , style=Style(color='k', donor='b-', acceptor='r-') )
      
      # all plots are lines
      plot( ... , ... , style='-' )
      
      # all plots are green
      plot( ... , ... , style=Style(color='g') )
      
      # color cycles for every line drawn (like matplotlib default)
      plot( ... , ... , style=Style(color='auto') )
  """
  ITEMS = ('donor', 'acceptor', 'fret', 'trap', 'total_counts')
  STYLES = (':','-','-',':','--')

  TRAP_STYLE_NAME = 'trap'

  def __init__(self, **styles):
    color = styles.pop('color', None)
    if color is None:
      color = next_color_in_cycle()
    elif color == 'auto':
      color = '' # lets matplotlib decide
    default_style = [color+linestyle for linestyle in Style.STYLES]
    for item,style in styles.items():
      if style and not _has_color(style):
        styles[item] = color + style
    super(Style, self).__init__(
      zip(Style.ITEMS, default_style),
      **styles)

  @classmethod
  def with_default_style(cls, style):
    style_dict = dict(zip(cls.ITEMS, [style]*len(cls.ITEMS)))
    return cls(**style_dict)
    
  def __setitem__(self, key, value):
    if key not in Style.ITEMS:
      raise ValueError('Item "{}" is not plotted and has no style')
    super(Style, self).__setitem__(key, value)


def save(fname=None, **kwargs):
    fname = fname or plt.gca().get_title()
    plt.savefig(fname, transparent=True, **kwargs)

def stack(p, style=None, **kwargs):
  trap = p.trap
  fret = p.fret
  plot(fret, pull=trap, style=style, **kwargs)

def stacks(pulls, style=None, **options):
  options.setdefault('legend', None)
  for p in pulls:
    # label = labels.pop(0) if len(labels)>0 else ''
    stack(p, style=style, hold=True, label=p.filename, **options)

def counts_from_fret(fretdata, style=Style(), rows=2):
  time, donor, acceptor, fret = fretdata
  plt.hold(True)
  plt.subplot(rows,1,1)
  dstyle = style['donor']
  plt.plot(time, donor, dstyle, label='Donor')
  plt.plot(time, acceptor, style['acceptor'], label='Acceptor')
  plt.plot(time, acceptor+donor, style['total_counts'], label='Total')
  plt.ylabel("Counts")
  plt.ylim(-100, max(acceptor+donor+100))
  plt.xlim(xmin=time[0],xmax=time[-1])
  plt.xticks([])

  plt.subplot(rows,1,2)
  plt.ylabel('FRET')
  plt.plot(time, fret, style['fret'])
  plt.xlabel('Time')
  plt.xlim(xmin=time[0],xmax=time[-1])
  plt.ylim(-0.05,1.05)
  return plt.gcf().get_axes()

def fret(p, style=Style(), rows=2):
  return counts_from_fret(p.fret, style=style, rows=rows)

def frets(exps, rows=2):
  plt.hold(True)
  for p in exps:
    fret(p, rows=rows)

def fec(p, style=':', *args, **options):
  options.setdefault('label', p.filename)
  out = plt.plot(p.trap.ext, p.trap.f, style, *args, **options)
  plt.xlabel('Extension (nm)')
  plt.ylabel('Force (pn)')
  plt.gca().autoscale_view(tight=True)
  return out

def fecs(exps, style=':', *args, **options):
  plt.hold(True)
  for p in exps:
    fec(p, style, *args, **options)

def setdefaults(dict_, **defaults):
  for name,default in defaults.iteritems():
    dict_.setdefault(name, default)
  return dict_

def subplots(nrows=1):
  return [plt.subplot(nrows,1,row) for row in range(1,nrows+1)]

def twinxaxis():
  return plt.gca(), plt.gca().twinx()

def fretfec(p, stylefec='-', stylefret='o:', axes=(), fec_opts={}, fret_opts={}):
  if any(axes):
    axfret,axforce = axes
  else:
    axfret, axforce = twinxaxis()
  fret,ext = to_fret_ext(p)
  plt.hold(True)
  fec_opts.setdefault('axes', axforce)
  setdefaults(fret_opts, **fretfec.fret_opts)
  setdefaults(fec_opts, **fretfec.fec_opts)
  fec(p, stylefec, **fec_opts)
  axfret.plot(ext, fret, stylefret, **fret_opts)
  plt.xlabel('Extension (nm)')
  plt.ylabel('FRET efficiency')
  axforce.autoscale_view(tight=True)
  return axfret, axforce

fretfec.fret_opts = {
  'markersize': 9
}
fretfec.fec_opts = {
  'linewidth': 2.5
}

from collections import deque

def segmented(p, stylefec='-', stylefret='o:', every=1, boundaries=[], colors=[], axes=(), exps=[], **options):
  axfret, axforce = axes if any(axes) else subplots(2)
  stylefret=stylefret or stylefec
  linewidth = options.pop('linewidth', 3)
  ratio = p['sampling_ratio']
  trap = p.trap
  fret, ext = to_fret_ext(p)
  plt.hold(True)
  colors = deque(colors)
  for x in exps:
    axforce.plot(x.trap.ext, x.trap.f, 'k:', lw=1.5)
  for start_n in range(len(fret)):
    if start_n == 0 or start_n in boundaries or (not boundaries and start_n%every==0):
      color = colors.popleft() if colors else next_color_in_cycle()
    fret_segment = fret[start_n:start_n+2]
    ext_segment = ext[start_n:start_n+2]
    axfret.plot(ext_segment, fret_segment, stylefret, color=color, markeredgewidth=0, **options)
    fecext,fecf = trap[start_n*ratio:(start_n+1)*ratio+1].fec
    axforce.plot(fecext, fecf, stylefec, color=color, markeredgewidth=0, lw=linewidth, **options)
  axfret.set_ylim(-0.05,1.05)
  axforce.autoscale_view(tight=True)
  axfret.set_xlim(axforce.get_xlim())
  axforce.set_ylabel('Force (pN)')
  axfret.set_ylabel('FRET efficiency')
  axfret.set_xlabel('Extension (nm)')

def _has_color(style_string):
  return style_string and len(style_string) > 0 and (style_string[0].isalnum() or style_string[0]=='#')


def plot(data, pull=None, style=None, **kwargs):
  """Plot FretData and/or TrapData as stacked subplots of counts, FRET, and FEC/FDC
  @data: datatypes.AbstractData
  @pull: datatypes.AbstractData
  @style: dict or str
  """
  if isinstance(style, str):
    style = Style.with_default_style(style)
  elif isinstance(style, dict):
    style = Style(**style)
  else:
    style = Style()
  loc = kwargs.pop('legend', 'best')
  title = kwargs.pop('title','')
  label = kwargs.pop('label', '')
  displayFRET = kwargs.pop('show_fret', True) and hasattr(data,'fret')
  kwargs.setdefault('markersize', 2)

  hold=kwargs.pop('hold', None)
  if hold is not None:
    plt.hold(hold)

  if not pull and hasTrapData(data):
    pull = TrapData.fromObject(data)

  num = 0
  if displayFRET:
    num += 2
  if pull:
    num += 1
    FEC = kwargs.pop('FEC', num==1)
  else:
    FEC = kwargs.pop('FEC', False)

  if num == 0:
    raise ValueError("Don't know how to plot arguments: need TrapData or FretData")

  if data is None and pull is not None:
    layout = iter([(num,1,num)])
  else:
    layout = iter((num,1,x) for x in range(1,num+1))

  ax1 = None
  donor, acceptor = 'donor', 'acceptor'
  if displayFRET:
    plt.subplot(*next(layout))
    not hold and plt.cla()
    plt.hold(True)
    ax1 = subplot(data.time, data.donor, style[donor], label=donor, **kwargs)
    subplot(data.time, data.acceptor, style[acceptor], label=acceptor, axes=('Time (s)','Counts'), **kwargs)
    plt.hold(hold)
    if loc is not None:
      plt.legend(loc=loc,ncol=2,prop={'size':'small'})
  if displayFRET:
    subplot(data.time, data.fret, style['fret'], layout=next(layout), 
              axes=('Time (s)','FRET'), **kwargs)
    plt.ylim(-0.1, 1.1)

  ax2 = None
  if pull:
    trap_style = Style.TRAP_STYLE_NAME
    x_coord,x_label = (pull.ext,'Extension (nm)') if FEC else (pull.sep,'Separation (nm)')
    ax2 = subplot(x_coord, pull.f, style[trap_style], layout=next(layout), axes=(x_label,'Force (pN)'), label=label, **kwargs)
    if loc is not None:
      plt.legend(loc=loc,ncol=2,prop={'size':'small'})

  first_plot = ax1 or ax2
  if title:
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

def hist(data, bins=50, hold=False, **kwargs):
  plt.hold(hold)
  counts,bins,patches = plt.hist(data, bins=bins, **kwargs)
  delta = bins[1]-bins[0]
  bins = (bins-delta/2)[1:]
  return bins,counts
