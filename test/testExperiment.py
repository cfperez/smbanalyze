import unittest
import operator
import os.path as path
import os
import pickle
from mock import Mock, patch, MagicMock
from numpy import array

from smbanalyze import experiment, fileIO, datatypes

class FigStub(object):
  def __init__(self):
    self.plotCalled = False
  def plot(self, *args, **kwargs):
    self.plotCalled = True
  def annotate(self, text, location):
    self.annotateCalled = True

def attrgetter(attr):
  return lambda x: map(operator.attrgetter(attr), x)

filegetter = lambda f: attrgetter('filename')(f)

def setUp():
  global pulls, LOADED_FILES
  os.chdir('test/')
  with patch('smbanalyze.experiment.Figure') as mock:
    pulls = experiment.fromMatch('test')
  LOADED_FILES = map( path.normpath, 
    [r'test_s1m1', r'test_s1m2', r'test_s1m3', ]
  )

def testRipFitting():
  pull = pulls[2]
  pull.fitHandles(800, 8)
  assert pull.fitRip(987)['Lc1'] > 10 #== 16.475625684642072

def testHandleFitting():
  pass

def testExtensionOffset():
  pass

def testPullingFromFile():
  pullLoad = [ experiment.Pulling.fromFile(f) for f in LOADED_FILES ]
  assert filegetter(pullLoad) == LOADED_FILES

def testPullingLoadimg():
  pull = pulls[1] # test_s1m2
  img = pull.loadimg()
  assert path.splitext(img.filename)[0] == path.splitext(pull.filename)[0]

def testExperimentFromMatch():
  filenames = filegetter(pulls)
  assert set(map(path.normpath, filenames)) == set(LOADED_FILES)

def testExperimentFromMatchWithOpenLoopExperiments():
  pass

def testExperimentFromFiles():
  loaded = experiment.fromFiles(LOADED_FILES)
  assert filegetter(loaded) == LOADED_FILES

def testPullingSave():
  pull = pulls[0]
  fname = 'save_s1m1.exp'
  with patch('smbanalyze.experiment.pickle.dump') as mock:
    pull.save(fname)
    mock.assert_called()

def testPullingLoad():
  fname = 'save_s1m1.exp'
  with patch('smbanalyze.experiment.pickle.load') as mock:
    mock.return_value = pulls[0]
    pull = experiment.Pulling.load(fname)
    mock.assert_called()
    assert type(pull) == experiment.Pulling
    assert pull.trap == pulls[0].trap
    assert pull.fret == pulls[0].fret

def testExperimentPlot():
  for a_pull in pulls:
    a_pull.figure = FigStub()
    a_pull.plot()
    assert a_pull.figure.plotCalled == True

def testList():
  pull_list = experiment.List(pulls)
  assert pull_list[0] == pulls[0]
  matched = pull_list.matching('s1m2')
  assert len(matched) == 1
  assert matched[0] == pulls[1]
  assert isinstance(matched, experiment.List)

def testListAggregatePull():
  pull_list = experiment.List(pulls)
  aggregated = pull_list.aggregate('trap')
  rows, cols = aggregated.shape
  total_rows = sum(p.shape[0] for p in pull_list.get('trap'))
  assert total_rows == rows

def testListAggregatePullNoFRET():
  pull_list = experiment.List(pulls).not_has('fret')
  aggregated = pull_list.aggregate('trap')
  assert type(aggregated) == datatypes.TrapData
  rows, cols = aggregated.shape
  total_rows = sum(p.shape[0] for p in pull_list.get('trap'))
  assert total_rows == rows
  
def testListAggregatePullEmpty():
  pass

def testListGet():
  pull_list = experiment.List(pulls)
  info = pull_list.get('filename')
  assert set(info) == set(LOADED_FILES)

def testListCall():
  pull_list = experiment.List(pulls)
  def testFunc():
    return True
  for p in pull_list:
    p.testFunc = testFunc
  assert pull_list.call('testFunc') == [True]*len(pull_list)

def testHasAnyAttrReturnsTrue():
  mock = Mock()
  assert experiment.hasAnyAttr(mock, 'test','foo','bar')

def testHasAnyAttrReturnsFalse():
  mock = Mock(spec=['foo'])
  assert experiment.hasAnyAttr(mock, 'test', 'bar') == False

def tearDown():
  pass

class TestDatatypes(unittest.TestCase):
  def _startPatch(self, to_patch):
    patch_obj = patch(to_patch)
    self.addCleanup(patch_obj.stop)
    return patch_obj.start()

  def setUp(self):
    def load_file_mock(filename, **kwargs):
      return {}, array([[1,2,3],[4,5,6]])

    self.load = self._startPatch('smbanalyze.datatypes.load')

  def testLoad(self):
    data = array([[1,2,3],[4,5,6]])
    self.load.return_value = ({},data)
    fdata = datatypes.FretData.fromFile('testing')
    self.load.assert_called_with('testing', comments=fileIO.toSettings)
    self.assertEqual(data.tolist(), fdata.data.tolist())

class TestLoadingExperimentsFromFile(unittest.TestCase):

  def _startPatch(self, to_patch):
    patch_obj = patch(to_patch)
    self.addCleanup(patch_obj.stop)
    return patch_obj.start()

  def setUp(self):
    self.exists = self._startPatch('smbanalyze.experiment.path.exists')
    self.exists.return_value = True
    self.fdata = self._startPatch('smbanalyze.experiment.FretData')
    data = [[1,2,3],[1,2,3],[1,2,3],[1,2,3]]
    self.fdata.fromFile.return_value = datatypes.FretData.fromFields(*data)

    self.tdata = self._startPatch('smbanalyze.experiment.TrapData')
    self.tdata.fromFile.return_value = datatypes.TrapData.fromFields(*data[:3])

  def testFromFileReturnsPullingExpWithStrFileNoTime(self):
    filename = 'construct_0.5nM_s1m1'
    self.checkExperimentWithFilename(filename, experiment.Pulling)
    self.tdata.fromFile.assert_called_with(filename+fileIO.PULL_FILE)
    self.fdata.fromFile.assert_called_with(filename+fileIO.FRET_FILE)

  def testFromFileReturnsOpenLoopExpWithForceInBasename(self):
    filename = 'construct_0.5nM_s1m1_5pN'
    self.checkExperimentWithFilename(filename, experiment.OpenLoop)
    self.tdata.fromFile.assert_called_with(filename+fileIO.PULL_FILE)
    self.fdata.fromFile.assert_called_with(filename+fileIO.FRET_FILE)

  def testFromFileReturnsOpenLoopExpWithTimeInBasename(self):
    filename = 'construct_0.5nM_s1m1_2min'
    self.checkExperimentWithFilename(filename, experiment.OpenLoop)
    self.tdata.fromFile.assert_called_with(filename+fileIO.PULL_FILE)
    self.fdata.fromFile.assert_called_with(filename+fileIO.FRET_FILE)

  def testFromFileReturnsOpenLoopExpWithStrFileAndTime(self):
    filename = 'construct_0.5nM_s1m1_2min.str'
    self.checkExperimentWithFilename(filename, experiment.OpenLoop)
    self.tdata.fromFile.assert_called_with(filename)

  def testFromFileReturnsOpenLoopExpWithFretFileAndTime(self):
    filename = 'construct_0.5nM_s1m1_2min.fret'
    self.checkExperimentWithFilename(filename, experiment.OpenLoop)
    self.fdata.fromFile.assert_called_with(filename)

  def testFromFileReturnsOpenLoopExpWithBasenameNoStrFile(self):
    filename = 'construct_0.5nM_s1m1_2min'
    self.exists.return_value = False
    self.tdata.fromFile.return_value = None
    self.checkExperimentWithFilename(filename, experiment.OpenLoop)

  def checkExperimentWithFilename(self, filename, exp_type):
    exp = experiment.fromFile(filename)
    self.assertIsInstance(exp, exp_type)
    self.assertEqual(exp.filename, fileIO.splitext(filename)[0])
    self.assertEqual(exp.trap, self.tdata.fromFile.return_value)
    self.assertEqual(exp.fret, self.fdata.fromFile.return_value)


class TestOpenLoopLoading(unittest.TestCase):

  def setUp(self):
    self.hasFiletype = experiment.OpenLoop.hasFiletype
    self.force_time_test = lambda f,t: 'SJF4_0.5nM_s1m3_3_{}pN_2{}.fret'.format(f,t)
    
    self.time_types = ('min','s')
    self.forces = (0,1,5.5,12)

  def testhasFiletypeReturnsTrueWithForceTimeInFilename(self):
    for force, ttype in zip(self.forces, self.time_types):
      self.checkForceTimeInFilename(force, ttype)

  def checkForceTimeInFilename(self, force, ttype):
    self.assertTrue( self.hasFiletype(self.force_time_test(force,ttype)) )

  def testhasFiletypeReturnsTrueWithForceInFilename(self):
    for force in self.forces:
      self.checkForceInFilename(force)

  def checkForceInFilename(self, force):
    force_test = 'SJF4_0.5nM_s1m3_3_{}pN.fret'.format(force)
    self.assertTrue( self.hasFiletype(force_test) )

  def testhasFiletypeReturnsTrueWithTimeInFilename(self):
    for tt in self.time_types:
      self.checkTimeTest(tt)

  def checkTimeTest(self, time_type):
    time_test = 'SJF4_0.5nM_s1m3_3_2{}.fret'.format(time_type)
    self.assertTrue( self.hasFiletype(time_test) )

  @patch('smbanalyze.experiment.path.exists')
  @patch('smbanalyze.experiment.TrapData')
  @patch('smbanalyze.experiment.FretData')
  def testOpenLoopFromFile(self, fdata, tdata, exists_mock):
      exists_mock.return_value = False
      tdata.fromFile.side_effect = IOError()
      data = [[1,2,3],[1,2,3],[1,2,3],[1,2,3]]
      data = datatypes.FretData.fromFields(*data)
      fdata.fromFile.return_value = data
      fname = 'testing'
      exp = experiment.OpenLoop.fromFile(fname)
      fdata.fromFile.assert_called_with(fname+fileIO.FRET_FILE)
      self.assertEquals(exp.fret, data)
