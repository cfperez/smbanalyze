'''Test smbanalyze.datatypes
'''
from nose.tools import *
from smbanalyze.datatypes import FretData, TrapData

def setup():
	global trap, trapmeta
	trapmeta = {'outer': {'inner':1}}
	trap = TrapData([[1,2,3],[2,3,4],[3,4,5],[4,5,6]], trapmeta)

def test_TrapData_column_from_attribute():
	for index,field in enumerate(TrapData._fields):
		assert_items_equal(getattr(trap,field), trap.data[:,index])

def test_TrapData_metadata_nested_access():
	eq_(trapmeta['outer']['inner'], trap.metadata['outer.inner'])
