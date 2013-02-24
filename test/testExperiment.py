import operator
from smbanalyze import Experiment, FileIO

def attrgetter(attr):
  return lambda x: map(operator.attrgetter(attr), x)

filegetter = lambda f: attrgetter('file')(f)

def setUp():
  global pulls, LOADED_FILES
  pulls = Experiment.fromMatch('test/test')
  LOADED_FILES = [r'test\test_s1m1', r'test\test_s1m2', r'test\test_s1m3']

def testRipFitting():
  pull = pulls[2]
  pull.fitHandles(800, 9.09)
  assert pull.fitRip(987) == 16.475625684642072

def testHandleFitting():
  pass

def testFECfitting():
  pass

def testPullingLoad():
  pullLoad = [
    Experiment.Pulling.fromFile(r'test\test_s1m1.str'),
    Experiment.Pulling.fromFile(r'test\test_s1m2'),
    Experiment.Pulling.fromFile(r'test\test_s1m3'),]
  assert map(repr, pullLoad) == map(repr, pulls)

def testPullingLoadWithFret():
  pass

def testExperimentFromMatch():
  filenames = filegetter(pulls)
  assert filenames == LOADED_FILES

def testExperimentFromFiles():
  loaded = Experiment.fromFiles(LOADED_FILES)
  assert filegetter(loaded) == LOADED_FILES

def testFromFile():
  pass

def tearDown():
  pass
