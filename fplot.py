import matplotlib.pyplot as plt

from types import PullFretData,FretData,hasPullData,hasFretData

def plotall(objList, **kwargs):
  for obj in objList:
    plt.figure()
    obj.plot(**kwargs)

def plot(data, pull=None, **kwargs):
  loc = kwargs.get('legend', 'best')
  title = kwargs.get('title','')
  FEC = kwargs.get('FEC',False)
  displayFRET = kwargs.get('show_fret', True)

  hold=kwargs.get('hold', None)
  if hold is not None:
    plt.hold(hold)

  if pull and not hasPullData(data):
    data = PullFretData(*(pull+data))

  num = kwargs.get('numplot',subplotsNeeded(data))
  if not displayFRET: num -= 1
  if num==0:
    raise ValueError("Don't know how to plot argument: maybe missing named fields")

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

  if hasPullData(data):
    x_coord,x_label = (data.ext,'Extension (nm)') if FEC else (data.sep,'Separation (nm)')
    _subplot(x_coord, data.f, '.', layout=next(layout), axes=(x_label,'Force (pN)'))

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
