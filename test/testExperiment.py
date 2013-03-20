import operator
import stubble
from smbanalyze import experiment, fileIO, fplot

class FigStub(object):
  __metaclass__ = stubble.stubclass(experiment.Figure)
  def __init__(self):
    self.plotCalled = False
  def plot(self, *args, **kwargs):
    self.plotCalled = True

def attrgetter(attr):
  return lambda x: map(operator.attrgetter(attr), x)

filegetter = lambda f: attrgetter('filename')(f)

def setUp():
  global pulls, LOADED_FILES
  pulls = experiment.fromMatch('test/test')
  LOADED_FILES = [r'test\test_s1m1', r'test\test_s1m2', r'test\test_s1m3']

def testRipFitting():
  pull = pulls[2]
  pull.fitHandles(800, 9.09)
  assert pull.fitRip(987)['Lc1'] == 16.475625684642072

def testHandleFitting():
  pass

def testFECfitting():
  pass

def testPullingLoad():
  pullLoad = [
    experiment.Pulling.fromFile(r'test\test_s1m1.str'),
    experiment.Pulling.fromFile(r'test\test_s1m2'),
    experiment.Pulling.fromFile(r'test\test_s1m3'),]
  assert map(repr, pullLoad) == map(repr, pulls)

def testPullingLoadWithFret():
  pass

def testExperimentFromMatch():
  filenames = filegetter(pulls)
  assert filenames == LOADED_FILES

def testExperimentFromFiles():
  loaded = experiment.fromFiles(LOADED_FILES)
  assert filegetter(loaded) == LOADED_FILES

def testExperimentPlot():
  for a_pull in pulls:
    a_pull.figure = FigStub()
    a_pull.plot()
    assert a_pull.figure.plotCalled == True

def testFromFile():
  pass

def tearDown():
  pass
