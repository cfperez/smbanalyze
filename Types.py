from numpy import all

import collections

TYPE_CHECKING='STRICT'

FretData_fields = ('time','donor','acceptor','fret')
PullData_fields = ('ext','f','sep')

FretData = collections.namedtuple('FretData', FretData_fields)
PullData = collections.namedtuple('PullData', PullData_fields)
PullFretData = collections.namedtuple('PullingFretData',PullData_fields+FretData_fields)

def _hasData(datatype):
  fields = getattr(datatype,'_fields',datatype)
  if TYPE_CHECKING=='STRICT':
    return lambda obj: all(map(hasattr,[obj]*len(fields),fields))
  else:
    return lambda obj: len(fields)==len(obj)

hasFretData = _hasData(FretData)
hasPullData = _hasData(PullData)
hasPullFretData = _hasData(PullFretData)
