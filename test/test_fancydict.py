from nose.tools import *
from smbanalyze.fancydict import nesteddict

def test_nesteddict_update_preserves_1level_nested_dicts():
	d=nesteddict({'a':1})
	d.update( {'b': {'test':'inner'}} )
	d.update(c={'test':'inner'})
	eq_(d['c.test'], 'inner')
	eq_(d['b.test'], 'inner')

def test_nesteddict_fromdict_creates_nesteddict():
	d = nesteddict.from_dict(
		{'a': 1,
		'b': {'test':'inner'}
		})
	eq_(d['a'], 1)
	eq_(d['b.test'], 'inner')
	eq_(d['b'], {'test': 'inner'})

def test_nesteddict_fromdict_creates_2level_nesteddict():
	d = nesteddict.from_dict(
		{'a': 1,
		'b': {'test':{'inner':1}}
		})
	eq_(d['a'], 1)
	assert_is_instance(d['b'], nesteddict)
	eq_(d['b.test.inner'], 1)

def test_nesteddict_init_preserves_1level_nested_dicts():
	d = nesteddict(
		{'a': 1,
		'b': {'test':'inner'}
		})
	eq_(d['a'], 1)
	eq_(d['b.test'], 'inner')
	eq_(d['b'], {'test': 'inner'})

def test_nesteddict_setitem_creates_nested_key():
	d = nesteddict(
		{'a': 1})
	d['test.test'] = True
	ok_(d['test.test'])

def test_nesteddict_from_dict_uses_dot_access():
	ok_(nesteddict.from_dict({'test.test': True})['test.test'])