import unittest
from mock import Mock, MagicMock, patch

from smbanalyze import fileIO as fio

glob_mock = Mock(spec='glob.glob')

@patch('glob.glob', glob_mock)
class TestFlist(unittest.TestCase):

  def testCallsWithGlobOperator(self):
    fio.flist('SJ')
    glob_mock.assert_called_with('*SJ*')
    
  def testMultipleGlobs(self):
    fio.flist('SJ','s1m1')
    glob_mock.assert_called_with('*SJ*s1m1*')

class TestParseFilename(unittest.TestCase):
  def setUp(self):
    self.template = '{construct}_0.5nM_s{slide}m{mol}_{pull}_{force}pN_{min}min.fret'

  def filenameFromParams(self, **params):
    return self.template.format(**params)

  def testReturnsFilenameData(self):
    params = dict(
      construct='SJF4', slide=1, mol=3,
      pull=3, force=0, min='2')

    finfo = fio.parseFilename(
      self.filenameFromParams(**params))

    self.checkReturnData(params, finfo)

  def testReturnsDataWithFloatingPointForce(self):
    params = dict(
      construct='SJF4', slide=1, mol=3,
      pull=3, force=10.5, min='2')
    finfo = fio.parseFilename(
      self.filenameFromParams(**params))
    self.checkReturnData(params, finfo)

  def checkReturnData(self, input, output):
    for param,value in input.items():
      self.assertEquals(str(getattr(output, param)), str(value))

  def testReturnsNoneWithBadFilename(self):
    fail = 'SJF4_0.5nM_s1m3_2hours.fret'
    finfo = fio.parseFilename(fail)
    self.assertEquals(finfo, None)
