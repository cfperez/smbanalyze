''' api.py
'''
from pymongo import MongoClient
from operator import itemgetter
from smbanalyze.experiment import List, Pulling
from bson.binary import Binary
from cPickle import dumps, loads

__all__ = ['lessthan', 'lessorequal']


def to_binary(nparray, subtype=128):
    return Binary(dumps(nparray, protocol=2), subtype)

def from_binary(binary):
    return loads(binary)

from pymongo.son_manipulator import SONManipulator
def fix_dict_keys(dict_):
    return {key.replace('.','_'): val for key,val in p.metadata.iteritems()}

class TransformToPickle(SONManipulator):
    def transform_incoming(self, son, collection):
        for (key,val) in son.items():
            if isinstance(val, Pulling):
                son[key] = [to_binary(val.trap), to_binary(val.fret), to_binary(val.metadata)]
            elif isinstance(val, dict):
                son[key] = self.transform_incoming(val, collection)
        return son
    
    def transform_outgoing(self, son, collection):
        for (key,val) in son.items():
            if isinstance(val,list) and getattr(val[0], 'subtype', 0) == 128:
                son[key] = Pulling(
                    trap=from_binary(val[0]),
                    fret=from_binary(val[1]),
                    **from_binary(val[2]))
            elif isinstance(val,dict):
                son[key] = self.transform_outgoing(val, collection)
        return son

CLIENT = None
def connect(database='data'):
    global CLIENT
    CLIENT = CLIENT or MongoClient()
    db = getattr(CLIENT, database)
    db.add_son_manipulator(TransformToPickle())
    return db

class FindResults(object):
    def __init__(self, cursor):
        self._list = list(cursor)

    def all(self):
        return self._list

    def select(self, *fields):
        return FindResults(itemgetter(*fields)(p) for p in self._list)
    def __iter__(self):
        return iter(self._list)

def select(cursor, *fields):
    return (itemgetter(*fields)(p) for p in cursor)

def find(db, **search):
    return db.find({key:val for key,val in search.iteritems()})

def get_exp(db, **search):
    return List(p['experiment'] for p in find(db, **search))

def exp_to_dict(exp):
    to_insert = {key.replace('.','_'): val for key,val in exp.metadata.iteritems()}
    to_insert['experiment'] = exp
    return to_insert

def save_exp(db, exp, extra={}, **extra_):
    to_insert = exp_to_dict(exp)
    to_insert.update(extra)
    to_insert.update(extra_)
    db.save(to_insert)
    
def save_all(db, exps):
    for p in exps:
        save_exp(db, p)
        
def lessthan(val):
    return {'$lt': val}
def lessorequal(val):
    return {'$lte': val}