from nose.tools import raises
from smbanalyze.image import ROI, ROIError

def setUp():
  global BL, TR, roi
  BL = (35, 50)
  TR = (80,100)
  roi = ROI( BL, TR, origin='absolute')

def testROItoRelative():
  old = ROI.copy(roi)
  origin = (30, 32)
  relROI = roi.toRelative( origin )
  assert str(old) == str(roi)
  assert relROI.origin == 'relative'
  assert old.left - relROI.left == origin[0]
  assert old.top - relROI.top == origin[1]

@raises(ROIError)
def testROItoRelativeFail():
  relROI = roi.toRelative( (100,100) )

@raises(ROIError)
def testROIFail():
  roi = ROI( (85,50), (80,100) )

def testROItoAbsolute():
  origin = (30, 32)
  absROI = roi.toAbsolute(origin)
  assert str(absROI) == str(roi)
