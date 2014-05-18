__all__ = ['upsample', 'downsample', 'down_average']

from numpy import sum

def upsample(arr, ratio):
    return [_ for r in arr for _ in range(r*ratio, (r+1)*ratio)]

def downsample(arr, ratio, offset=None):
    if offset is None:
        offset = ratio-1
    elif offset >= ratio:
        raise ValueError(
        '"offset" argument (%d) must be less than downsample "ratio"' % offset)
    size = len(arr)
    return arr[offset:size-(size%ratio):ratio]

def down_average(arr, ratio):
  return sum((downsample(arr, ratio, offset=i) for i in range(ratio)), axis=0)/ratio
