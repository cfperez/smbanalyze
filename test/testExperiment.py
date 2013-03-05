import operator
from smbanalyze import experiment, fileIO

def attrgetter(attr):
  return lambda x: map(operator.attrgetter(attr), x)

filegetter = lambda f: attrgetter('filename')(f)

def setUp():
  global pulls, LOADED_FILES
  pulls = experiment.fromMatch('test/test')
  LOADED_FILES = [r'test/test_s1m1', r'test/test_s1m2', r'test/test_s1m3']

def testRipFitting():
  pull = pulls[2]
  pull.fitHandles(800, 9.09)
  assert pull.fitRip(987) == 16.475625684642072

def testHandleFitting():
  pass

def testFECfitting():
  pass

def testPullingLoad():
  pullLoad = [ experiment.Pulling.fromFile(f) for f in LOADED_FILES ]
  assert filegetter(pullLoad) == LOADED_FILES

def testPullingLoadWithFret():
  pass

def testexperimentFromMatch():
  filenames = filegetter(pulls)
  assert set(filenames) == set(LOADED_FILES)

def testexperimentFromFiles():
  loaded = experiment.fromFiles(LOADED_FILES)
  assert filegetter(loaded) == LOADED_FILES

def testFromFile():
  pass

def tearDown():
  pass
