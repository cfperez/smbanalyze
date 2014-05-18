from numpy import where, iterable, diff, array
from matplotlib.mlab import find
from operator import isCallable
from datatypes import FretData
from collections import namedtuple
from fancydict import dotdict

ON_STATE = 1
OFF_STATE = 0


def find_states(arr, threshold):
    "Return array of binary states using simple thresholding"
    return where(arr > threshold, ON_STATE, OFF_STATE)


def count_states(state_array, state):
    "Return the number of states of given value in state_array"
    assert iterable(state_array)
    return len(find(state_array == state))


def find_transitions(state_array):
    "Return array of transitions: 0 no transition, 1 pos transition, -1 neg transition"
    assert iterable(state_array)
    deriv_array = diff(state_array)
    return deriv_array


def count_transitions(transition_array, transition_direction):
    "Return the number of transitions with the given direction"
    assert transition_direction == ON_STATE - \
        OFF_STATE or transition_direction == OFF_STATE - ON_STATE
    return count_states(transition_array, transition_direction)


def state_lifetime(state_array, state):
    "Return array of lifetimes for specified state"
    assert iterable(state_array)
    count = 0
    lifetime = []
    for s in state_array:
        if s == state:
            count += 1
        elif count > 0:
            lifetime += [count]
            count = 0
    if count > 0:
        lifetime += [count]
    return array(lifetime)


def count_blinking(exp, acc_threshold, donor_threshold):
    "Returns # of acceptor, donor blinks as found by simple threshold"
    transitions = lambda data, threshold: find_transitions(
        find_states(data, threshold))
    acc = count_transitions(transitions(exp.fret.acceptor, acc_threshold), -1)
    donor = count_transitions(transitions(exp.fret.donor, donor_threshold), 1)
    return acc, donor

Lifetime = namedtuple('Lifetime', 'on off')

def single(threshold, filter_=None):
    return StateFinder(threshold, filter_)

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

    def __init__(self, threshold, filter_=None):
        if filter_ is not None or isCallable(filter_):
            raise ValueError('"filter_" arg must be callable')
        self.threshold = threshold
        self.filter = filter_

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

    def lifetimes(self, input_data): # state_array=None):
        state_array = self.findStates(input_data)
        on = state_lifetime(state_array, ON_STATE)
        off = state_lifetime(state_array, OFF_STATE)
        return Lifetime(on, off)

    def findTransitions(self):
        return find_transitions(self.states)


class StateFinderMultipleThreshold(StateFinder):

    def __init__(self, *thresholds):
        self.thresholds = thresholds
        self.filter = None

    @classmethod
    def with_filter(cls, filter_, *thresholds):
        filtered = cls(*thresholds)
        filtered.filter = filter_
        return filtered

    def findStates(self, input_data, thresholds=[]):
        thresholds = thresholds or self.thresholds
        return tuple(find_states(input_data, threshold) for threshold in thresholds)

        for attr in self.threshold:
            self.states[attr] = find_states(
                getattr(data, attr), self.threshold[attr])
        return self.states

    def lifetimes(self, state_dict=None, attr=None):
        if state_dict:
            assert isinstance(state_dict, dict)
        if attr:
            raise ValueError('WARNING: Use of an attribute name in function call has been DEPRECATED.\n'
                             + 'Call lifetimes() with no arguments.')
        state_dict = state_dict or self.states
        return dotdict([(attr_in_dict, super(StateFinderMultipleThreshold, self).lifetimes(states))
                        for (attr_in_dict, states) in state_dict.iteritems()])

    def findTransitions(self):
        self.transitions = {}
        for key in self.states:
            self.transitions[key] = find_transitions(self.states[key])
        return self.transitions


class StateFinderFretData(StateFinderMultipleThreshold):

    def __init__(self, **thresholds): #donor=[], acceptor=[], fret=[]):
        self.thresholds = {attr:val if iterable(val) else [val] for attr,val in thresholds.items()} #dict(donor=donor, acceptor=acceptor, fret=fret)

    def findStates(self, fret_data):
        self._meta_data = fret_data.metadata
        self._frame_time = fret_data.metadata['exposurems'] / 1000.
        return {attr: super(StateFinderFretData, self).findStates(getattr(fret_data,attr), thresholds) 
            for attr,thresholds in self.thresholds.items()}
        return super(StateFinderFretData, self).findStates(fret_data)

    def lifetimes(self, fret_data): #state_dict=None, attr=None):
        # assert isinstance(state_dict, (dict, type(None)))
        # t = self._frame_time
        states = self.findStates(fret_data)
        # lifetimes = super(StateFinder, self).lifetimes() for state_array in 
        lifetime = super(StateFinderFretData, self).lifetimes(state_dict, attr)
        for key, val in lifetime.iteritems():
            lifetime[key] = Lifetime(*(v * t for v in val))
        return lifetime
