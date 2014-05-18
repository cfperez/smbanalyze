''' api.py
'''
from pymongo import MongoClient
from operator import itemgetter
from smbanalyze.experiment import ExpList, Pulling
from smbanalyze.datatypes import TrapData, FretData
from bson.binary import Binary
from cPickle import dumps, loads
from smbanalyze import date

__all__ = ['_id', 'set_', 'lessthan', 'lessorequal', 'get_exp',
 'connect', 'select', 'find', 'save_exp', 'save_all', 'copy', 'update']

CLIENT = None
def connect(database='data'):
    global CLIENT
    CLIENT = CLIENT or MongoClient()
    db = getattr(CLIENT, database)
    #db.add_son_manipulator(TransformToPickle())
    return db

def copy(from_, to_, **search):
    '''Insert every item in collection from_ into to_'''
    for item in from_.find(search):
        to_.insert(item)

def data_to_db(p, dtype):
    data = getattr(p, dtype)
    d = data.metadata.copy()
    d['data'] = data.data.tolist()
    return d

def exp_to_db(p):
    to_save = p.metadata.copy()
    if p.trap:
        # to_save['trap']['data'] = p.trap.data.tolist()
        to_save['trap'] = data_to_db(p, 'trap')
    if p.fret:
        to_save['fret'] = data_to_db(p, 'fret')
    return to_save

def db_to_exp(p):
    trap = p['trap'].pop('data')
    trap_data = TrapData(trap, **p['trap'])
    if 'fret' in p:
        fret = p['fret'].pop('data', None)
        fret_data = FretData(fret, **p['fret']) if fret else None
    return Pulling(trap_data,fret_data, p)

def _id(d):
    if isinstance(d, str):
        return {'_id':d}
    else:
        try:
            return {'_id':d['_id']}
        except KeyError:
            raise ValueError(
                'Object %s has no key "_id" for database')

def to_binary(nparray, subtype=128):
    return Binary(dumps(nparray, protocol=2), subtype)

def from_binary(binary):
    return loads(binary)

from pymongo.son_manipulator import SONManipulator

def fix_dict_keys(dict_):
    return {key.replace('.','_'): val for key,val in dict_.iteritems()}

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


def closedb(db):
    db.connection.close()

class FindResults(object):
    def __init__(self, cursor):
        self._list = list(cursor)

    def count(self):
        return len(self._list)

    def all(self):
        return self._list

    def has(self, **conditions):
        # (p for p in self if )
        pass

    def select(self, *fields):
        def getter(p):
            return tuple(p.get(field,None) for field in fields)
        return FindResults(map(getter,self._list))

    def __iter__(self):
        return iter(self._list)

def select(cursor, *fields):
    return (itemgetter(*fields)(p) for p in cursor)

def find(db, **search):
    date_ = search.get('date', None)
    if date_ and isinstance(date_, tuple):
        search['date'] = date.date(*date_)
    return FindResults( db.find(search) )

def find_one(db, **search):
    return db.find_one(**search)

def get_exp(db, **search):
    return ExpList(db_to_exp(p) for p in find(db, **search))

def exp_to_dict(exp):
    to_insert = {key.replace('.','_'): val for key,val in exp.metadata.iteritems()}
    to_insert['experiment'] = exp
    return to_insert

def save_exp(db, exp, extra={}, **extra_):
    exp.metadata.update(extra, **extra_)
    to_insert = exp_to_db(exp)
    # to_insert.update(extra, **extra_)
    db.save(to_insert)
    
def save_all(db, exps, extra={}, **extra_):
    for p in exps:
        save_exp(db, p, extra, **extra_)
        
def update(db, obj, **fields):
    return db.update(_id(obj), set_(**fields))

def update_all(db, match, **fields):
    return db.update(match, set_(**fields), multi=True)

def set_(**fields):
    return {'$set': fields}

def in_(vals):
    return {'$in': vals}
    
def greaterthan(val):
    return {'$gt': val}
    
def lessthan(val):
    return {'$lt': val}

def lessorequal(val):
    return {'$lte': val}