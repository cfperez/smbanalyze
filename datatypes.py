import logging
from numpy import all, asarray, sign, vstack
from operator import isSequenceType, attrgetter
from collections import Iterable
from fileIO import load, toSettings, fileIOError
import constants
import copy

logger = logging.getLogger(__name__)
logger.setLevel(constants.logLevel)
logger.addHandler(constants.logHandler)

TYPE_CHECKING = 'STRICT'

FretData_fields = ('time', 'donor', 'acceptor', 'fret')
TrapData_fields = ('ext', 'f', 'sep')


class AbstractData(object):

    def __init__(self, data, **meta):
        data = asarray(data)
        try:
            if data is not None and data.shape[1] != len(self._fields):
                raise ValueError('Data argument must have fields for {}'.format(self._fields))
        except IndexError:
            raise ValueError('Data argument must be a two-dimensional array')
        self.data = data
        self._original_data = None
        self.metadata = {}
        self.metadata.update(meta)

    @classmethod
    def fromObject(cls, obj):
        try:
            return cls(copy.copy(obj.data), **obj.metadata)
        except AttributeError:
            raise ValueError(
                'Constructor method only takes object with AbstractData interface')

    @classmethod
    def name(cls):
        return cls.__name__.lower()

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
        me = cls(data, **meta)
        return me

    @classmethod
    def aggregate(cls, abstractdata, sort_by=None):
        assert isinstance(abstractdata, Iterable)
        assert len(abstractdata) > 0
        assert all(isinstance(d, AbstractData)
                   or d is None for d in abstractdata)
        if sort_by and sort_by not in cls._fields:
            raise ValueError(
                'sort_by argument must be a field in this data type')
        key = lambda e: e[cls._fields.index(sort_by)] if sort_by else None
        data = map(attrgetter('data'), abstractdata)
        return cls(asarray(sorted(vstack(data), key=key)))

    @classmethod
    def fromFields(cls, *args, **meta):
        if len(args) != len(cls._fields):
            raise ValueError(
                "Number of arguments to {0} must match _fields".format(cls.__name__))
        return cls(asarray(args).T, **meta)

    def copy(self):
        return self.__class__.fromObject(self)

    @property
    def shape(self):
        return self.data.shape

    def __len__(self):
        return len(self.data)

    def __eq__(self, other):
        return all(self.data == other.data) and self.metadata == other.metadata

    def __ne__(self, other):
        return not self == other

    METADATA_CHECK = {'trap': ('step_size', 'sampling_time'),
                      'fret': ('exposurems', 'frames', 'gain', 'binning')}

    def __add__(self, other):
        # if self.meta != other.meta:
            # logger.warning('M')
        return type(self)(vstack((self.data, other.data)))

    def at(self, **kwargs):
        if len(kwargs) > 1:
            raise ValueError('Use only one keyword representing a field')

        field, value = kwargs.popitem()
        if field not in self._fields:
            raise ValueError('Keyword argument must be a field')

        return self[getattr(self, field) > value][0]
        index = search_monotonic(getattr(self, field), value)
        return self[index]

    def __getattr__(self, attr):
        if attr in self._fields:
            attr_position = self._fields.index(attr)
            if len(self.shape) == 1:
                return self.data[attr_position]
            else:
                return self.data[:, attr_position]
        else:
            try:
                return getattr(super(AbstractData, self), attr)
            except AttributeError:
                raise AttributeError("{0} has no attribute '{1}'".format(
                    self.__class__.__name__, attr))

    @property
    def T(self):
        return self.data.T

    def __iter__(self):
        return iter(self.T)

    def __getitem__(self, key):
        return type(self)(self.data[key].view(), **self.metadata)

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


class TrapData(AbstractData):
    _fields = TrapData_fields

    def maskFromLimits(self, x=None, f=None, limits=()):
        if x is None and f is None:
            raise ValueError('Must specific either x limits of f limits')
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

    def select(self, x=None, f=None, limits=(0, -1)):
        return self[self.maskFromLimits(x, f, limits)]

    @property
    def fec(self):
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
        self.metadata['stiffness'] = new_k
        beadRadii = self.metadata.get('bead_radii', constants.sumOfBeadRadii)

        displacement = self.trap.f / min(current_k)
        ratio = 1 + min(new_k) / max(new_k)

        self.trap.f = displacement * min(new_k)
        self.trap.ext = self.trap.sep - beadRadii - \
            displacement * ratio - self._ext_offset
        return self


class FretData(AbstractData):
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
