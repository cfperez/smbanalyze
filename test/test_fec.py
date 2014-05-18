'''
Tests for fec.py
'''
import unittest
import mock
from smbanalyze.curvefit import Fit
from collections import OrderedDict

from smbanalyze import fec
from numpy import sqrt

def quadrature_error(*vals):
    return sqrt(sum(map(lambda x:x**2, vals)))

class BaseRipErrors(unittest.TestCase):
	def setUp(self):
		self.mock_fit = mock.Mock(spec=Fit)

	def set(self, error):
		self.mock_fit.error = OrderedDict(error)

	def check(self):
		return fec.rip_errors(self.mock_fit)

class TestRipErrors(BaseRipErrors):
	def test_returns_single_error_one_rip(self):
		errors_dict = OrderedDict(
			(('Lc',1.5), ('Lc1', 2.5))
		)
		self.set(errors_dict)
		errors = self.check()
		self.assertItemsEqual(errors, [errors_dict['Lc1']])

	def test_returns_errors_two_rips(self):
		self.set(
			(('Lc',1.1), ('Lc1',2.5), ('Lc2',3.7))
		)
		out = self.check()
		print out
		calculated = [2.5,quadrature_error(3.7,2.5)]
		print calculated
		self.assertItemsEqual(out, calculated)