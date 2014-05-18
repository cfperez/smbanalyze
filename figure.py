import matplotlib.pyplot as plt
from fileIO import DEFAULT_FIGURE_EXT

class Figure(object):
  def __init__(self, fig_id=None):
    self.figure_id = fig_id
    self._figure = None
    self.reversed = False

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
    self.makeCurrent()
    #plt.figure(self.figure_id)
    return self.figure

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
    return self._figure

  def pickPoints(self, num_of_pts=2):
  	return self.figure.ginput(num_of_pts)

  def pickRegions(self, num=1):
    points = sorted(x for x,f in self.pickPoints(num*2))
    return [(points[i],points[i+1]) for i in range(0,len(points),2)]   
    
  def makeCurrent(self):
    if not self.exists:
      raise RuntimeError('Figure is not visible')
    plt.figure(self.figure_id)
    return self

  def xlim(self, xmin=None, xmax=None, reverse=False):
    assert not reverse or not (xmin or xmax)
    if reverse and not self.reversed:
      xmax, xmin = plt.xlim()
      self.reversed = True
    return plt.xlim(xmin, xmax)

  def plot(self, *args, **kwargs):
    self.show()
    plt.plot(*args, **kwargs)
  #   try:
  #     # Treat the first argument as an object that can plot itself...
  #     return args[0].plot(*args[1:], **kwargs)
  #   except AttributeError,ValueError:
  #     # ...unless it can't
  #     return plot(*args, **kwargs)

  # def plotall(self, *args, **kwargs):
  #   self.show()
  #   plotall(*args, **kwargs)

  def clear(self):
    if self.visible:
      self._figure.clf()
      self.reverse = False
      plt.draw()
    return self
      
  def close(self):
    if self.exists:
      plt.close(self.figure_id)
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
        filename += DEFAULT_FIGURE_EXT
    else:
      filename = 'Figure {0}{1}'.format(self.figure.number, DEFAULT_FIGURE_EXT)
    self._figure.set_size_inches(*size)
    self._figure.savefig(filename, bbox_inches='tight', pad_inches=0.1)
