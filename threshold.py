from numpy import where, iterable, diff, array
from matplotlib.mlab import find

ON_STATE = 1
OFF_STATE = 0

def find_states(arr, threshold):
  return where(arr>threshold, ON_STATE, OFF_STATE)

def count_states(state_array, state):
  assert iterable(state_array)
  return len( find(state_array==state) )

def find_transitions(state_array):
  assert iterable(state_array)
  deriv_array = diff(state_array)
  return deriv_array

def count_transitions(transition_array, transition_direction):
  assert transition == ON_STATE-OFF_STATE or transition == OFF_STATE-ON_STATE
  return count_states(transition_array, transition_direction)

def state_lifetime(state_array, state):
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
  sf = StateFinder.fromThreshold(donor=1500, acceptor=1200, fret=0.5)
  sf.findStates( long_trace.fret )
  
  sf.states['donor']
  sf.states['acceptor']
  sf.states['fret']
  """

  @classmethod
  def fromThreshold(cls, **thresholds):
    states = cls()
    states.thresholds = thresholds
    return states

  def __call__(self, data):
    return self.findStates(data)

  def findStates(self, data):
    self.states = {}
    for attr in self.thresholds:
      self.states[attr] = find_states( getattr(data, attr), self.thresholds[attr] )
    return self.states

  def lifetime(self, attr):
    if attr not in self.states:
      raise ValueError('"{}" not in states'.format(attr))
    state_array = self.states[attr]
    on = state_lifetimes(state_array, ON_STATE)
    off = state_lifetimes(state_array, OFF_STATE)
    return on, off

  def findTransitions(self):
    self.transitions = {}
    for key in self.states:
      self.transitions[key] = find_transitions(self.states[key])
    return self.transitions