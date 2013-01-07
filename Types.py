import collections

FretData = collections.namedtuple('FretData',('time','donor','acceptor','fret'))
PullData = collections.namedtuple('PullData', ('ext','f','sep'))
PullFretData = collections.namedtuple('PullingFretData',PullData._fields+FretData._fields)
