import unittest
from nose.tools import raises
from smbanalyze.image import ROI, ROIError

class ROItests(unittest.TestCase):
  def setUp(self):
    self.BL = (35, 50)
    self.TR = (80,100)
    self.roi = ROI( self.BL, self.TR, origin='absolute')

  def testInit(self):
    self.assertEqual(self.BL, (self.roi.left,self.roi.bottom))
    self.assertEqual(self.TR, (self.roi.right,self.roi.top))

  def testROItoRelative(self):
    old = ROI.copy(self.roi)
    origin = (30, 32)
    relROI = self.roi.toRelative( origin )
    assert relROI.origin == 'relative'
    assert old.left - relROI.left == origin[0]
    assert old.top - relROI.top == origin[1]

  @raises(ROIError)
  def testROItoRelativeFail(self):
    relROI = self.roi.toRelative( (100,100) )

  @raises(ROIError)
  def testROIFail(self):
    roi = ROI( (85,50), (80,100) )

  def testROItoAbsolute(self):
    origin = (30, 32)
    absROI = self.roi.toAbsolute(origin)
    assert str(absROI) == str(self.roi)
