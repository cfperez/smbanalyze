import operator
import matplotlib.pyplot as plt
from numpy import concatenate

from datatypes import TrapData,hasTrapData,hasFretData

def plotall(objList, **kwargs):
  labels = kwargs.pop('labels', [])
  for obj in objList:
    label = labels.pop(0) if len(labels)>0 else ''
    plot(obj, label=label, **kwargs)
  legend = kwargs.get('legend', 2)
  plt.legend(loc=legend)

def plot(data, pull=None, **kwargs):
  loc = kwargs.get('legend', 'best')
  title = kwargs.get('title','')
  label = kwargs.get('label', '')
  FEC = kwargs.get('FEC',False)
  displayFRET = kwargs.get('show_fret', hasattr(data,'fret'))

  hold=kwargs.pop('hold', None)
  if hold is not None:
    plt.hold(hold)

  if not pull and hasTrapData(data):
    pull = TrapData.fromObject(data)

  num = 0
  if hasattr(data, 'donor'): num += 1
  if displayFRET and hasFretData(data): num += 1
  if pull: 
    num += 1
    if num == 1:
      FEC = kwargs.setdefault('FEC', True)

  if num == 0:
    raise ValueError("Don't know how to plot arguments: need TrapData or FretData")
  #if not displayFRET and data is not None: num = 1
  #elif hasFretData(data): num = 2
  #else: num = 0 #raise ValueError("Don't understand plot options")
  #if pull: 
  #  num += 1
  #if num==0:
  #  raise ValueError("Don't know how to plot argument: missing named fields")

  layout = iter((num,1,x) for x in range(1,num+1))

  ax1 = None
  if hasFretData(data):
    plt.subplot(*next(layout))
    not hold and plt.cla()
    plt.hold(True)
    ax1 = _subplot(data.time, data.donor, label='donor')
    _subplot(data.time, data.acceptor, label='acceptor', axes=('','counts'))
    plt.hold(hold)
    plt.legend(loc=loc,ncol=2,prop={'size':'small'})
    if displayFRET:
      _subplot(data.time, data.fret, layout=next(layout), 
                axes=('Seconds','FRET'))

  ax2 = None
  if pull:
    x_coord,x_label = (pull.ext,'Extension (nm)') if FEC else (pull.sep,'Separation (nm)')
    ax2 = _subplot(x_coord, pull.f, '.', layout=next(layout), axes=(x_label,'Force (pN)'), label=label)

  first_plot = ax1 or ax2
  first_plot.set_title(title)
  plt.show()

def subplotsNeeded(data):
  num = 0
  if hasattr(data, 'donor') and hasattr(data, 'acceptor'):
    num += 1
  if hasattr(data,'fret'):
    num += 1
  if hasPullData(data):
    num += 1
  return num

def _subplot(*args,**kwargs):
  sub = kwargs.pop('layout',())
  axes = kwargs.pop('axes', ())
  hold = kwargs.get('hold', None)

  if sub:
    ax = plt.subplot(*sub)
  else:
    ax = plt.gca()
  if hold is not None:
    plt.hold(hold)
  plt.plot(*args,**kwargs)
  ax.autoscale_view(tight=True)
  if axes:
    try:
      plt.xlabel(axes[0])
      plt.ylabel(axes[1])
    except IndexError:
      raise ValueError('_subplot expects labels for BOTH axes')
  return ax

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
