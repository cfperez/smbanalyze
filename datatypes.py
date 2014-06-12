import logging
from numpy import all, asarray, sign, vstack, ndarray, arange
from operator import and_
from operator import isSequenceType, attrgetter
from collections import Iterable
from fileIO import load, toSettings, fileIOError
import constants
import copy
from fancydict import nesteddict
from contextlib import contextmanager

logger = logging.getLogger(__name__)
logger.setLevel(constants.logLevel)
logger.addHandler(constants.logHandler)

TYPE_CHECKING = 'STRICT'

FretData_fields = ('time', 'donor', 'acceptor', 'fret')
TrapData_fields = ('ext', 'f', 'sep')

class Mask(ndarray):
    '''Not used! Here as an example of subclassing ndarray'''
    def __new__(cls, bool_array, *args):
        assert isinstance(bool_array, ndarray)
        obj = asarray(bool_array, dtype='bool').view(cls)
        return obj

    def above(self, above_func):
        pass

@contextmanager
def transform_error(original_error, to_error):
    try:
        yield
    except original_error:
        raise to_error

_index_to_column_error = lambda i: transform_error(IndexError, ValueError("Data does not have column %d" % i))

def _field_properties(fields):
    '''Create name,property tuples to access column i in data (DataType)
    '''
    # Must define fget/fset in this way OUTSIDE of for loop
    def property_funcs(i):
        def fget(self):
            with _index_to_column_error(i):
                return self.data[i].view() if len(self.data.shape)==1 else self.data[:,i].view()
        def fset(self, value):
            with _index_to_column_error(i):
                self.data[:,i] = value
        return fget,fset

    for index,name in enumerate(fields):
        fget, fset = property_funcs(index)
        yield name, property(fget, fset, doc="Field property %s (%d)" % (name,index))

class DataType(type):
    """Metaclass for building datatypes with _fields attribute access"""
    def __new__(cls, name, bases, clsdict):
        # Add the properties (e.g. 'ext','f','sep' for TrapData) to the class before creating it
        clsdict.update(
            _field_properties(clsdict.get('_fields',[]))
        )
        obj = super(DataType, cls).__new__(cls, name, bases, clsdict)
        return obj


class Data(object):
    __metaclass__ = DataType

    def __init__(self, data, meta={}):
        data = asarray(data)
        if data is not None and not self._is_data_shape_ok(data.shape):
            logger.warning('TrapData should have fields for {}'.format(self._fields))
        self.data = data
        self._original_data = None
        self.metadata = nesteddict.from_dict(meta)

    @classmethod
    def _is_data_shape_ok(self, shape):
        field_size = len(self._fields)
        return len(shape) == 1 and shape[0] == field_size or shape[1] == field_size

    @classmethod
    def fromObject(cls, obj):
        try:
            return cls(copy.copy(obj.data), obj.metadata)
        except AttributeError:
            raise ValueError(
                'Constructor method only takes object with Data interface')

    @classmethod
    def fromFile(cls, filename):
        try:
            meta, data = load(filename, comments=toSettings)
        except fileIOError as e:
            if not e.isError:
                print e.strerror
                data = load(filename)
                meta = {}
            else:
                raise
        else:
            meta['filename'] = filename
        me = cls(data, meta)
        return me

    @classmethod
    def aggregate(cls, dataiter, sort_by=None):
        assert isinstance(dataiter, Iterable)
        assert len(dataiter) > 0
        assert all(isinstance(d, Data)
                   or d is None for d in dataiter)
        if sort_by and sort_by not in cls._fields:
            raise ValueError(
                'sort_by argument must be a field in this data type')
        key = lambda e: e[cls._fields.index(sort_by)] if sort_by else None
        data = map(attrgetter('data'), dataiter)
        return cls(asarray(sorted(vstack(data), key=key)))

    flatten = aggregate

    @classmethod
    def fromFields(cls, *args, **meta):
        if len(args) != len(cls._fields):
            raise ValueError(
                "Number of arguments to {0} must match _fields".format(cls.__name__))
        return cls(asarray(args).T, meta)

    def copy(self):
        return self.__class__.fromObject(self)

    @property
    def shape(self):
        return self.data.shape

    def __len__(self):
        return len(self.data)

    def __eq__(self, other):
        if hasattr(other, 'data'):
            return all(self.data == other.data) and self.metadata == other.metadata
        else:
            return self.data == other

    def __ne__(self, other):
        return not self == other

    def __add__(self, other):
        return type(self)(vstack((self.data, other.data)))

    def where(self, *conditions):
        '''Returns data where ALL conditions are true
        '''
        return self[reduce(and_, conditions)]

    def at(self, **kwargs):
        return self.where(*[getattr(self,field)>=value for field,value in kwargs.iteritems()])[0]

    @property
    def T(self):
        return self.data.T

    def __iter__(self):
        return iter(self.T)

    def __getitem__(self, key):
        if len(self.data.shape) == 1:
            return self.data[key]
        return type(self)(self.data[key].view(), self.metadata)

    def __repr__(self):
        return repr(self.data)

    @classmethod
    def _normalizeLimits(cls, limits, min_max, assume_max_limit=True):
        if limits is None:
            return min_max
        elif not isSequenceType(limits):
            limits = [limits]
        if len(limits) == 2:
            return limits
        elif len(limits) == 1:
            if assume_max_limit:
                return min_max[0], limits[0]
            else:
                return limits[0], min_max[1]


class TrapData(Data):
    _fields = TrapData_fields

    @property
    def time(self):
        if 'sampling_time' not in self.metadata:
            raise AttributeError(
                'sampling_time must be set in metadata to calculate time')
        return arange(1,len(self))*self.metadata['sampling_time']
        
    def maskFromLimits(self, x=None, f=None, limits=()):
        if x is None and f is None:
            raise ValueError('Must specify either x limits or f limits')
        if limits:
            start, stop = limits
            ext_fit, f_fit = self.ext[start:stop], self.f[start:stop]
        else:
            ext_fit, f_fit = self.ext, self.f

        min_f, max_f = TrapData._normalizeLimits(f,
                                                 min_max=(
                                                     min(f_fit), max(f_fit)),
                                                 assume_max_limit=True
                                                 )

        min_ext, max_ext = TrapData._normalizeLimits(x,
                                                     min_max=(
                                                         min(ext_fit), max(ext_fit)),
                                                     assume_max_limit=False
                                                     )

        between = lambda s, a, b: (s >= a) & (s <= b)
        return between(ext_fit, min_ext, max_ext) & between(f_fit, min_f, max_f)

    def mask_from_interval(self, ext, f=None):
        return self.maskFromLimits(ext, f)

    def mask_above(self, above):
        '''Return 2 masks: True if above above(ext)'''
        assert callable(above)
        above = self.f > above(self.ext)
        return above

    def make_masks(self, intervals):
        return map(self.mask_from_interval, intervals)

    def select(self, x=None, f=None, ext=None, limits=(0, -1)):
        x = x or ext
        return self[self.maskFromLimits(x, f, limits)]

    @property
    def fec(self):
        """Return (ext,force) data
        """
        x, f, s = self
        return x, f

    def meanStiffness(self):
        inverseAverage = lambda args: 1 / sum(map(lambda x: 1. / x, args))
        return inverseAverage(self.metadata.get('stiffness', constants.stiffness))

    def adjustOffset(self, ext=None, force=None):
        if ext:
            self.ext -= ext
        if force:
            self.force += force

    def recalculate(self, stiffness=None):
        if stiffness and len(stiffness) != 2:
            raise ValueError('Stiffness must be 2-tuple')
        current_k = self.metadata.get('stiffness', stiffness)
        new_k = stiffness or current_k
        self.metadata['stiffness'] = tuple(new_k)
        beadRadii = self.metadata.get('bead_radii', constants.sumOfBeadRadii)

        displacement = self.f / min(current_k)
        ratio = 1 + min(new_k) / max(new_k)

        self.f = displacement * min(new_k)
        self.ext = self.sep - beadRadii - displacement * ratio
        return self


class FretData(Data):
    _fields = FretData_fields

    def maskFromLimits(self, time, limits=(0, -1)):
        return

    def select(self, time=None):
        return self[self.maskFromLimits(time)]

def search_monotonic(ar, value):
    shifted = ar - value
    if value <= ar[0]:
        return 0
    elif value >= ar[-1]:
        return -1

    start_sign = sign(shifted[0])
    for n, (current, last) in enumerate(and_prev(shifted)):
        if sign(current) != start_sign:
            return n if abs(current) < abs(last) else min(n - 1, 0)
    return -1

def and_prev(iterable, default=None):
    last = default
    for x in iterable:
        yield x, last
        last = x


def _hasData(datatype):
    fields = getattr(datatype, '_fields', datatype)
    if TYPE_CHECKING == 'STRICT':
        return lambda obj: all(map(hasattr, [obj] * len(fields), fields))
    else:
        return lambda obj: len(fields) == len(obj)

hasFretData = _hasData(FretData)
hasTrapData = _hasData(TrapData)
hasTrapFretData = lambda x: hasFretData(x) and hasTrapData(x)
