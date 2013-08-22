from smbanalyze.threshold import StateFinder, StateFinderMultipleThreshold
from smbanalyze.datatypes import FretData
import unittest
from numpy import array

class TestStateFinder(unittest.TestCase):

    def setUp(self):
        self.threshold = 1000
        self.level_change_data = array([0, 0, 1, 1, 0])
        self.input_data = self.level_change_data*self.threshold*2
        self.SF = StateFinder.fromThreshold(self.threshold)
        self.output = self.SF(self.input_data)

    def test_fromThreshold_returns_StateFinder(self):
        self.assertIsInstance(self.SF, StateFinder)

    def test_findStates_returns_correct_states(self):
        data = self.level_change_data * self.threshold * 2
        self.assertListEqual(
            self.SF(data).tolist(), self.level_change_data.tolist())


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


class TestStateFinderFretData(TestStateFinderMultipleThreshold):
    pass