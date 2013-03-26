import matplotlib.pyplot as plt

from datatypes import TrapData,hasTrapData,hasFretData

def plotall(objList, **kwargs):
  for obj in objList:
    obj.plot(**kwargs)

def plot(data, pull=None, **kwargs):
  loc = kwargs.get('legend', 'best')
  title = kwargs.get('title','')
  FEC = kwargs.get('FEC',False)
  displayFRET = kwargs.get('show_fret', hasattr(data,'fret'))

  hold=kwargs.get('hold', None)
  if hold is not None:
    plt.hold(hold)

  if hasattr(data, '_to_plot'):
    args, kwargs = data._to_plot()
    return _subplot(*args, **kwargs)

  if not pull and hasTrapData(data):
    pull = TrapData.fromObject(data)

  if not displayFRET and data: num = 1
  elif hasFretData(data): num = 2
  else: num = 0
  if pull: num += 1
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
    plt.legend(loc=loc,ncol=2,prop={'size':'small'})
    if displayFRET:
      _subplot(data.time, data.fret, layout=next(layout), 
                axes=('Seconds','FRET'))

  if pull:
    x_coord,x_label = (pull.ext,'Extension (nm)') if FEC else (pull.sep,'Separation (nm)')
    _subplot(x_coord, pull.f, '.', layout=next(layout), axes=(x_label,'Force (pN)'))

  plt.gcf().get_axes()[0].set_title(title)
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
