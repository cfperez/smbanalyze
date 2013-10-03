import unittest
from nose.tools import raises
import operator
import os.path as path
import os
from mock import Mock, patch, MagicMock
from numpy import array, ndarray

from smbanalyze import experiment, fileIO, datatypes, fplot, curvefit

class TestCase(unittest.TestCase):
  def _startPatch(self, to_patch, **kwargs):
    patch_obj = patch(to_patch, **kwargs)
    self.addCleanup(patch_obj.stop)
    return patch_obj.start()

def attrgetter(attr):
  return lambda x: map(operator.attrgetter(attr), x)

filegetter = lambda f: attrgetter('filename')(f)

def setUp():
  global pulls, LOADED_FILES
  os.chdir('test/')
  experiment.Options['loading']['filename_matching'] = False
  with patch('smbanalyze.fplot.Figure'):
    pulls = experiment.fromMatch('test')
  LOADED_FILES = map( path.normpath, 
    [r'test_cond_s1m1', r'test_cond_s1m2', r'test_cond_s1m3']
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

@patch('smbanalyze.image.fromFile')
def testPullingLoadimg(fromFile):
  pull = pulls[1] # test_s1m2
  pull.loadimg()
  fromFile.assert_called_with(fileIO.add_img_ext(pull.filename))
  
@patch('smbanalyze.image.fromFile')
@raises(experiment.ExperimentError)
def testPullingLoadimgNoImage(fromFile):
  pull = pulls[1]
  fromFile.side_effect = IOError
  pull.loadimg()

def testExperimentFromMatch():
  filenames = filegetter(pulls)
  assert set(map(path.normpath, filenames)) == set(LOADED_FILES)

def testExperimentFromFiles():
  loaded = experiment.fromFiles(LOADED_FILES)
  assert filegetter(loaded) == LOADED_FILES

def testPullingSave():
  pull = pulls[0]
  fname = 'save_s1m1.exp'
  with patch('smbanalyze.experiment.pickle.dump') as mock:
    pull.save(fname)
    mock.assert_called()

@patch('smbanalyze.experiment.pickle.dump')
class ExperimentSave(object):
    """save() in Pulling and List"""

    def __init__(self, *args, **kwargs):
        self.fname = 'test'
        super(ExperimentSave, self).__init__(*args, **kwargs)

    def testPullingSaveFilenameEndsWithExp(self, dump):
        self.exp.save(self.fname)
        (p,fh),kw = dump.call_args
        self.assertEqual(fh.name, self.fname+'.exp')
        self.assertEqual(p, self.exp)

    def testPullingSaveEmptyFilenameEndsWithExp(self, dump):
        self.exp.metadata['filename'] = self.fname
        self.exp.save()
        self.assertEqual(self.fname+'.exp', self._filename_from_dump(dump))

    def _filename_from_dump(self, dump):
        (p,fh),kw = dump.call_args
        return fh.name

class TestPullingSave(ExperimentSave, unittest.TestCase):
    def setUp(self):
        self.exp = experiment.Pulling(None, None)

class TestOpenLoopSave(ExperimentSave, unittest.TestCase):
    def setUp(self):
        self.exp = experiment.OpenLoop(None, None)


def testPullingLoad():
  fname = 'save_s1m1.exp'
  with patch('smbanalyze.experiment.pickle.load', autospec=True) as mock:
    mock.return_value = pulls[0]
    pull = experiment.Pulling.load(fname)
    assert mock.called
    assert type(pull) == experiment.Pulling
    assert pull.trap == pulls[0].trap
    assert pull.fret == pulls[0].fret

def testExperimentPlot():
  for a_pull in pulls:
    a_pull._figure = Mock(autospec=fplot.Figure)
    a_pull.plot()

def testList():
  pulls = [Mock(autospec=experiment.Pulling)]*3
  pull_list = experiment.List(pulls)
  assert pull_list == pulls
  #matched = pull_list.matching('s1m2')
  #assert len(matched) == 1
  #assert matched[0] == pulls[1]
  #assert isinstance(matched, experiment.List)

def testListAggregatePull():
  pull_list = experiment.List(pulls)
  aggregated = datatypes.TrapData.aggregate(
    pull_list.get('trap'))
  rows, cols = aggregated.shape
  total_rows = sum(p.shape[0] for p in pull_list.get('trap'))
  assert total_rows == rows

def testListAggregatePullNoFRET():
  pull_list = experiment.List(pulls).not_has('fret')
  aggregated = datatypes.TrapData.aggregate(
    pull_list.get('trap'))
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

def tearDown():
  pass

class TestFretData(unittest.TestCase):
  def _startPatch(self, to_patch):
    patch_obj = patch(to_patch)
    self.addCleanup(patch_obj.stop)
    return patch_obj.start()

  def setUp(self):
    self.data = array([[1,2,3,4],[4,5,6,7]])
    self.load = self._startPatch('smbanalyze.datatypes.load')
    self.load.return_value = ({},self.data)

  def test_1D_data_attr_access(self):
    data = self.data[0]
    fdata = datatypes.FretData(data)
    self.assertEqual(fdata.shape, (4,))
    self.assertItemsEqual(data, fdata)
    for i,field in enumerate(fdata._fields):
      self.assertEqual(data[i], getattr(fdata, field))

  def testLoad(self):
    fdata = datatypes.FretData.fromFile('test_cond_s1m1')
    self.load.assert_called_with('test_cond_s1m1', comments=fileIO.toSettings)
    self.assertEqual(self.data.tolist(), fdata.data.tolist())

  def test_getitem_returns_one_row_FretData(self):
    fdata = datatypes.FretData(self.data)
    for row in range(len(fdata)):
      first_row = fdata[row]
      self.assertIsInstance(first_row, datatypes.FretData)
      self.assertItemsEqual(first_row, self.data[row])

  def test_getitem_slice_returns_same_type(self):
    fdata = datatypes.FretData(self.data)
    self.assertIsInstance(fdata[:], datatypes.FretData)
    self.assertEqual(fdata, fdata[:], msg="Array slice doesn't match original")


class TestLoadingExperimentsFromFile(TestCase):

  def setUp(self):
    self.exists = self._startPatch('smbanalyze.experiment.opath.exists')
    self.exists.return_value = True
    self.fdata = self._startPatch('smbanalyze.experiment.FretData')
    data = [[1,2,3],[1,2,3],[1,2,3],[1,2,3]]
    self.fdata.fromFile.return_value = datatypes.FretData.fromFields(*data)

    self.tdata = self._startPatch('smbanalyze.experiment.TrapData')
    self.tdata.fromFile.return_value = datatypes.TrapData.fromFields(*data[:3])

  def checkFromFileCalls(self, filename):
    self.tdata.fromFile.assert_called_with(filename+fileIO.PULL_FILE)
    self.fdata.fromFile.assert_called_with(filename+fileIO.FRET_FILE)
    
  def testFromFileReturnsPullingExpWithStrFileNoTime(self):
    filename = 'construct_0.5nM_s1m1'
    self.checkExperimentWithFilename(filename, experiment.Pulling)
    self.checkFromFileCalls(filename)
    
  def testFromFileReturnsPullingExpWithStrFileNoConditions(self):
    filename = 'construct_100pM_s1m1'
    self.checkExperimentWithFilename(filename, experiment.Pulling)
    self.checkFromFileCalls(filename)

  def testFromFileReturnsOpenLoopExpWithForceInBasename(self):
    filename = 'construct_0.5nM_s1m1_5pN'
    self.checkExperimentWithFilename(filename, experiment.OpenLoop)
    self.checkFromFileCalls(filename)

  def testFromFileReturnsOpenLoopExpWithTimeInBasename(self):
    filename = 'construct_0.5nM_s1m1_2min'
    self.checkExperimentWithFilename(filename, experiment.OpenLoop)
    self.checkFromFileCalls(filename)

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
    self.hasFiletype = experiment.OpenLoop.filenameMatchesType
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

  def testhasFilteypeReturnsFalseWithPullingSyntax(self):
      self.assertFalse( self.hasFiletype('construct_100nM_s1m1'))
      self.assertFalse( self.hasFiletype('construct_100nM_s1m1_2'))
      
  def checkTimeTest(self, time_type):
    time_test = 'SJF4_0.5nM_s1m3_3_2{}.fret'.format(time_type)
    self.assertTrue( self.hasFiletype(time_test) )

  @patch('smbanalyze.experiment.opath.exists')
  @patch('smbanalyze.experiment.TrapData')
  @patch('smbanalyze.experiment.FretData')
  def testOpenLoopFromFile(self, fdata, tdata, exists_mock):
      exists_mock.return_value = False
      tdata.fromFile.side_effect = IOError()
      data = [[1,2,3],[1,2,3],[1,2,3],[1,2,3]]
      data = datatypes.FretData.fromFields(*data)
      fdata.fromFile.return_value = data
      fname = 'test_cond_s1m1'
      exp = experiment.OpenLoop.fromFile(fname)
      fdata.fromFile.assert_called_with(fname+fileIO.FRET_FILE)
      self.assertEquals(exp.fret, data)

class TestExperimentFromMatch(TestCase):
    def setUp(self):
        self.flist = self._startPatch('smbanalyze.fileIO.flist')
        self.openloop = self._startPatch('smbanalyze.experiment.OpenLoop.fromFile')
        self.pulling = self._startPatch('smbanalyze.experiment.Pulling.fromFile')
        experiment.Options.filtering.required_pulling_force = False
        experiment.Options.loading.filename_matching = True
        
    def testReturnsPullingAndOpenLoopExperiment(self):
        files = ['construct_100pM_s1m1_2','construct_10nM_s1m1_3_5pN']
        self.flist.return_value = map(fileIO.add_pull_ext, files)
        pulls = experiment.fromMatch('construct')
        self.assertEqual(len(pulls), 1)
        self.flist.assert_called_with('construct')
        self.pulling.assert_called_with(files[0])

class TestExperimentFindRip(TestCase):
  def setUp(self):
    tdata = datatypes.TrapData([[1000,5,1000],[1010,1,1010], [1020,2, 1020]])
    self.exp = experiment.Pulling(trap=tdata)

  @raises(experiment.ExperimentError)
  def test_error_on_no_trap_data(self):
    experiment.Pulling(None, None).findRip(900)

  @raises(experiment.ExperimentError)
  def test_error_on_no_fit(self):
    self.exp.findRip(900)

  def test_returns_1D_array_of_TrapData_fields(self):
    fit = MagicMock(autospec=curvefit.FitRegions)
    self.exp.fit = fit

class ListTest(unittest.TestCase):

    DATA_LENGTH = 10
    FILENAME = 'construct_100pM_s1m1_2'
 
    def setUp(self):
        self.trap_data = datatypes.TrapData.fromFields(
            [1]*ListTest.DATA_LENGTH, 
            [2]*ListTest.DATA_LENGTH,
            [3]*ListTest.DATA_LENGTH
            )
        self.to_collapse = [self.trap_data]*2

    def tearDown(self):
        pass

    def test_collapse_on_two_exp_with_trap_data(self):
        test_list = experiment.List(
            map(experiment.Pulling, self.to_collapse))
        collapsed = test_list.collapse()
        print collapsed.trap
        correct = datatypes.TrapData.aggregate(self.to_collapse)
        print correct
        self.assertItemsEqual(collapsed.trap.data.flat, 
            correct.data.flat)
        self.assertEqual(collapsed.filename, 'collapsed')

    def test_collapse_two_exp_with_filename(self):
        filename = 'test'
        test_list = experiment.List(
            map(experiment.Pulling, self.to_collapse))
        for p in test_list:
            p.metadata['filename'] = filename
        collapsed = test_list.collapse()
        self.assertEqual(collapsed.filename, filename+'_collapsed')         
