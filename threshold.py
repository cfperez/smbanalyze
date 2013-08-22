from numpy import where, iterable, diff, array
from matplotlib.mlab import find
from operator import isCallable
from datatypes import FretData
from collections import namedtuple
from smbanalyze.dotdict import dotdict

import unittest

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

    def lifetimes(self, state_array=None):
        if state_array is None:
            state_array = self.states
        on = state_lifetime(state_array, ON_STATE)
        off = state_lifetime(state_array, OFF_STATE)
        return Lifetime(on, off)

    def findTransitions(self):
        return find_transitions(self.states)


class TestStateFinder(unittest.TestCase):

    def setUp(self):
        self.threshold = 1000
        self.level_change_data = array([0, 0, 1, 1, 0])
        self.SF = StateFinder.fromThreshold(self.threshold)

    def test_fromThreshold_returns_StateFinder(self):
        self.assertIsInstance(self.SF, StateFinder)

    def test_findStates_returns_correct_states(self):
        data = self.level_change_data * self.threshold * 2
        self.assertListEqual(
            self.SF(data).tolist(), self.level_change_data.tolist())

    def test_lifetimes_returns_tuple(self):
        pass


class StateFinderMultipleThreshold(StateFinder):

    def __init__(self, **thresholds):
        self.threshold = thresholds

    def findStates(self, data):
        self.states = dotdict()
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


class TestStateFinderMultipleThreshold(unittest.TestCase):

    def setUp(self):
        self.SF = StateFinderMultipleThreshold(donor=1000)
        self.array = array([0, 0, 1, 1, 0])
        self.input_data = FretData.fromFields(*[self.array] * 4)

    def test_lifetimes_returns_dict_of_tuples(self):
        self.SF(self.input_data)
        lifetimes = self.SF.lifetimes()
        self.assertIsInstance(lifetimes, dict)
        for key in lifetimes:
            self.assertIsInstance(lifetimes[key], tuple)


class StateFinderFretData(StateFinderMultipleThreshold):

    def findStates(self, fret_data):
        self._meta_data = fret_data.metadata
        self._frame_time = fret_data.metadata['exposurems'] / 1000.
        return super(StateFinderFretData, self).findStates(fret_data)

    def lifetimes(self, state_dict=None, attr=None):
        assert isinstance(state_dict, (dict, type(None)))
        t = self._frame_time
        lifetime = super(StateFinderFretData, self).lifetimes(state_dict, attr)
        for key, val in lifetime.iteritems():
            lifetime[key] = Lifetime(*(v * t for v in val))
        return lifetime


class TestStateFinderFretData(TestStateFinderMultipleThreshold):
    pass
