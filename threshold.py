from numpy import where, iterable, diff, array
from matplotlib.mlab import find
from operator import isCallable
from datatypes import FretData

ON_STATE = 1
OFF_STATE = 0

def find_states(arr, threshold):
  "Return array of binary states using simple thresholding"
  return where(arr>threshold, ON_STATE, OFF_STATE)

def count_states(state_array, state):
  "Return the number of states of given value in state_array"
  assert iterable(state_array)
  return len( find(state_array==state) )

def find_transitions(state_array):
  "Return array of transitions: 0 no transition, 1 pos transition, -1 neg transition"
  assert iterable(state_array)
  deriv_array = diff(state_array)
  return deriv_array

def count_transitions(transition_array, transition_direction):
  "Return the number of transitions with the given direction"
  assert transition == ON_STATE-OFF_STATE or transition == OFF_STATE-ON_STATE
  return count_states(transition_array, transition_direction)

def state_lifetime(state_array, state):
  "Return array of lifetimes for specified state"
  assert iterable(state_array)
  count = 0
  lifetime = []
  for s in state_array:
    if s==state: count+=1
    elif count>0:
      lifetime += [count]
      count=0
  if count>0:
    lifetime += [count]
  return array(lifetime)

def count_blinking(exp, acc_threshold, donor_threshold):
  "Returns # of acceptor, donor blinks as found by simple threshold"
  transitions = lambda data, threshold: find_transitions(find_states(data, threshold))
  acc = count_transitions(transitions(exp.fret.acceptor, acc_threshold), -1)
  donor = count_transitions(transitions(exp.fret.donor, donor_threshold), 1)
  return acc, donor

class StateFinder(object):
  """
  Object for assigning states to single molecule trajectories.

  Create a StateFinder by calling the appropriate function then apply
  to the data of interest. It will calculate lifetimes and rates.

  sf = StateFinder.fromThreshold(donor=1500, acceptor=1200, fret=0.5)
  sf( p.fret )

  sf.states['donor']
  sf.states['acceptor']
  sf.states['fret']
  
  # OR 
  sf = StateFinder.fromThreshold(1300)
  sf( any_1D_array )
  """
  def __init__(self, threshold):
    self.threshold = threshold
    self.filter = None

  @classmethod
  def fromThreshold(cls, threshold=None, **thresholds):
    if threshold is not None:
      return StateFinder(threshold)
    else:
      return StateFinderMultipleThreshold(**thresholds)

  @classmethod
  def forFretData(cls, **thresholds):
    return StateFinderFretData(**thresholds)

    @classmethod
 
  def fromThresholdFiltered(cls, filter_func, threshold=None, **thresholds):
    assert isCallable(filter_func)
    raise NotImplementedError

  def __call__(self, data):
    return self.findStates(data)

  def addFilter(self, filter_func):
    assert isCallable(filter_func)
    self.filter = filter_func
  
  def findStates(self, input_data):
    "Calculate and return state array using threshold on data"
    data = self.filter(input_data) if self.filter else input_data
    self.states = find_states(data, self.threshold)
    return self.states

  def lifetime(self, state_array=None):
    if state_array is None:
      state_array = self.states
    on = state_lifetime(state_array, ON_STATE)
    off = state_lifetime(state_array, OFF_STATE)
    return on, off

  def rate(self, transition_array=None):
    pass

  def findTransitions(self):
    return find_transitions(self.states)

class StateFinderMultipleThreshold(StateFinder):
  def __init__(self, **thresholds):
    self.threshold = thresholds

  def findStates(self, data):
    self.states = {}
    for attr in self.threshold:
      self.states[attr] = find_states( getattr(data, attr), self.threshold[attr] )
    return self.states

  def lifetime(self, attr):
    if attr not in self.states:
      raise ValueError('"{}" not in states'.format(attr))
    return super(StateFinderMultipleThreshold, self).lifetime(self.states[attr])

  def findTransitions(self):
    self.transitions = {}
    for key in self.states:
      self.transitions[key] = find_transitions(self.states[key])
    return self.transitions

class StateFinderFretData(StateFinderMultipleThreshold):
  def findStates(self, fret_data):
    self._frame_time = fret_data.metadata['exposurems']/1000.
    return super(StateFinderFretData, self).findStates(fret_data)

  def lifetime(self, attr):
    t = self._frame_time
    on, off = super(StateFinderFretData, self).lifetime(attr)
    return on*t, off*t