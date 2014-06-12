'''
Test cases for experiment.py behaviors

- file loading
'''

from nose.tools import *
from smbanalyze.experiment import Pulling, to_fret_ext, ExpList
from smbanalyze.datatypes import FretData, TrapData
from numpy import all

def setup():
	global meta, trap, fret
	meta = dict(field1=1, field2='test')
	trap = TrapData([[1,2,3],[2,3,4],[3,4,5],[4,5,6]], meta)
	fret = FretData([[1,2,3,4],[2,3,4,5]], meta)

def test_pulling_trapdata():
	p = Pulling(trap)
	eq_(p.trap, trap)
	ok_(all(p.trap.data==trap.data))
	eq_(p.metadata['trap'], meta)
	# These fields have not been set and 
	# shouldn't be specified
	assert_not_in('fret.exposurems', p.metadata)
	assert_not_in('sampling_ratio', p.metadata)
	assert_is(p.fret, None)

def test_pulling_fretdata():
	p = Pulling(trap, fret)
	eq_(p.fret, fret)
	ok_(all(p.fret.data==fret.data))
	eq_(p.metadata['fret'], meta)

def test_pulling_to_fret_ext():
	p = Pulling(trap, fret)
	f,ext=to_fret_ext(p, 2)
	eq_(len(f), len(ext))

def test_experiment_list_has_fret_returns_list_with_fret():
	p = Pulling(trap, fret)
	mol = ExpList([p,p])
	has_fret = mol.has_attr('fret')
	eq_(has_fret, mol)
	print type(has_fret)
	assert_is_instance(has_fret, ExpList)

def test_experiment_list_nofret():
	p = Pulling(trap)
	mol = ExpList([p,p])
	has_fret = mol.has_attr('fret')
	eq_(has_fret, [])

def test_explist_collapses_trap_data():
	p = Pulling(trap, None, meta)
	mol = ExpList([p,p])
	collapsed = mol.collapse().trap
	aggregated = TrapData.aggregate([trap,trap], 'ext')
	print "Collapsed using ExpList.collapse():", collapsed.data
	print "Aggregated using TrapData.aggregate():", aggregated.data
	ok_(all(collapsed.data==aggregated.data), "Data is stacked and sorted by ext")
	eq_(aggregated.metadata, {}, "Metadata is ignored when aggregated")
