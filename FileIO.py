from __future__ import with_statement
import re
import collections
from datetime import datetime
import os
import operator
import glob
import ast
from inspect import isfunction

import numpy as np

from Types import FretData, PullData
from useful import toNum, toInt, isInt

IMAGE_FILE = '.img'
CAMERA_FILE = '.cam'
FRET_FILE = '.fret'
PULL_FILE = '.str'

REGISTERED_EXT = (IMAGE_FILE,CAMERA_FILE,FRET_FILE,PULL_FILE)

def flist(*globs):
  globs = list(globs)
  last = globs[-1]
  if isInt(last):
    globs[-1] = '_'+last

  return glob.glob('*%s*' % '*'.join(globs))

def load(fname, comments=False, header=False, **kwargs):
  base,extension = os.path.splitext(fname)
  fromLoad =  LOAD_FUNCTIONS[extension](fname, **kwargs)
  if isinstance(fromLoad, tuple):
    loadComments, loadHeader, loadData = fromLoad
    loadComments = comments(loadComments) if isfunction(comments) \
                    else loadComments
    output = [loadData]
    if header:
      output.insert(0,loadHeader)
    if comments:
      output.insert(0,loadComments)
    return tuple(output) if len(output)>1 else output[0]
  else:
    return fromLoad

def loadstr(fname, **kwargs):
  filecomments,fileheader,data = loaddat(fname,comments='#',**kwargs)
  if fileheader != ['extension','force','trapdistance']:
    raise IOError, "Stretch file must contain extension, force, and separation"
  return filecomments, fileheader, PullData(*data.T)
loadstr.extension=PULL_FILE

def toSettings(comments):
  # Careful! Using literal_eval to cast so that lines can contain
  # lists e.g. stiffness = [1,1.5]
  return parseSettingsFromLines(comments.splitlines(), ast.literal_eval)

def loaddat(filename, **kwargs):

  comment_line = kwargs.pop('comments', COMMENT_LINE)

  colnames = None
  comments = ''
  with open(filename,'rU') as fh:
    position = 0
    # loop stops after grabbing column names
    for number,line in enumerate(fh):
      if np.any( map(line.startswith, comment_line) ):
        comments += line.strip(' '+''.join(comment_line))
      elif line.isspace():
        continue
      elif colnames:
        position = number
        break
      elif any( map(str.isalnum, line.split()) ):
        colnames = line.lower().split()

    fh.seek(0)
    data = np.loadtxt(fh, skiprows=position, **kwargs)

  return comments, colnames, data

def loadfret(filename, **kwargs):
  comments,header,data = loaddat(filename,**kwargs)
  return comments, header, FretData(*data.T)
loadfret.extension=FRET_FILE

def savefret(filename, data, metadata=None, comments=''):
  if comments:
    comments += '\n'
  if metadata:
    comments += fromSettings(metadata)
  try:
    savedat(filename, (data.time,data.donor,data.acceptor,data.fret),
      header='time donor acceptor FRET', 
      fmt=('%.3f','%u','%u','%.5f'),
      comments=comments)
  except AttributeError:
    raise ValueError('Expects data with time, donor, acceptor, and fret attributes')

def loadsettings(filename, **kwargs):
  "Return dictionary of key/value pairs from LabView-style configuration file"
  cast = kwargs.get('cast') or str
  with open(filename) as f:
    return parseSettingsFromLines(f,cast)

def parseSettingsFromLines(lines, cast):
  heading = None
  settings = {}
  for line in strip_blank_and_comments(lines):
    header,key,value = parseLineToSetting(line)
    if not header:
      if heading: settings[heading][key] = cast(value)
      else: settings[key]=cast(value)
    else:
      heading=header
      settings[heading]={}
  return settings

def parseLineToSetting(line):
  m = re.match(r'\[(\w+)\]',line)
  if m:
    header = m.group(1).lower()
    return header,None,None
  else:
    key,value = line.split('=')
    return None,key.strip().lower().replace(' ','_'),value.strip()

def strip_blank_and_comments(iter_str):
  for line in strip_blank(iter_str):
    if not line.startswith('#'): yield line

def strip_blank(iter_str):
  for line in iter_str:
    line=line.strip()
    if line != '': yield line

def savedat(filename, data, header='', comments='', fmt='%.9e', delimiter='\t'):

  newline = '\n' if comments else ''
  header = ''.join(map(lambda line: COMMENT_LINE+' '+line, comments.splitlines(True))) \
      + newline + header.replace(' ',delimiter)
  if type(data) == tuple:
    data = np.array(data).T

  with open(filename, 'w') as fh:
    fh.write(header + '\n')
    if hasattr(data,'dtype') and data.dtype.names:
        fh.write( delimiter.join(data.dtype.names) + '\n' )
    np.savetxt(fh, data, fmt=fmt, delimiter=delimiter)

def loadimg(filename, datatype='>h', **kwargs):
  data = np.fromfile(filename, datatype)
  img_size = data[:3]
  img = np.delete( data, np.r_[:3, np.prod(img_size)+3:data.size] )
  try:
    img = img.reshape(img_size)
  except ValueError:
    raise IOError("Image file %s is corrupted, expected frame: %d, height: %d, width: %d" % 
        (filename, img_size[0], img_size[1], img_size[2]))
  return img

def loadcam(filename):
  "Return dictionary of camera settings as contained in .cam file"
  # Can be replaced by loadsettings and a handle for the DateTime key value
  # once the old-style .cam setting is supplanted

  settings = {}
  with open(filename, 'r') as f:
    for line in f.readlines():
      key,value = line.strip().split('\t')
      m = re.match('^\d+/', key)
      if key != 'DateTime' and not m:
        value = int(value)
      else:
        if m: # Handle old style cam files which don't have a DateTime key
          value = key + '\t' + value
          key = 'DateTime'
        value = datetime.strptime(value,'%m/%d/%Y %I:%M %p')

      settings[key.lower()] = value

  return settings
loadcam.extension=CAMERA_FILE

def savesettings(filename, file_mode, **settings):
  "save key/value object into LabView-style configuration file"

  with open(filename, file_mode) as f:
    f.write(fromSettings(settings))

## !! Skipping datetime key when saving because parseSettingFromLine()
## !! doesn't handle it yet
def fromSettings(settings):
  output = []
  HEADING_FMT = '[%s]\n'
  SETTING_FMT = '%s=%s\n'
  for header,subHead in settings.iteritems():
    try:
      settingsInHeading = subHead.iteritems()
      output += HEADING_FMT % header
      for setting,value in settingsInHeading:
        # UPDATE ME!
        if setting.lower() == 'datetime': continue
        output += SETTING_FMT % (setting,value)
    except AttributeError:
      output += SETTING_FMT % (header, subHead)
  return ''.join(output)

##################################################
## Filename Parsing
##################################################
def splitext(fname):
  basename,ext=os.path.splitext(fname)
  if ext not in REGISTERED_EXT:
    basename,ext=fname,''
  return basename,ext

change_extension = lambda x,y: os.path.splitext(x)[0]+y
add_img_ext = lambda x: x+IMAGE_FILE if not x.endswith(IMAGE_FILE) else x
add_cam_ext = lambda x: x+CAMERA_FILE if not x.endswith(CAMERA_FILE) else x
add_fret_ext = lambda x: x+FRET_FILE if not x.endswith(FRET_FILE) else x
add_pull_ext = lambda x: x+PULL_FILE if not x.endswith(PULL_FILE) else x

_named_ = lambda n,p: r'(?P<%s>%s)' % (n,p)
_opt_unit_ = lambda n,p,u: r'(?:_%s%s)?' % (_named_(n,p),u)
_opt_ = lambda n,p: _opt_unit_(n,p,'')
_opt_time_ = r'(?:_(?:(\d+)min)?(?:(\d+)s)?)?'
_build_pattern_ = lambda *x: ''.join(x)+'(?:\.\w+)?$'

_ALPHANUM_ = r'[a-zA-Z0-9]+'

FILENAME_SYNTAX = _build_pattern_( 
    '_'.join( [
      _named_('construct',_ALPHANUM_),
      _named_('conditions','.*'),
      r's(?P<slide>\d+)(?:m(?P<mol>\d+))?'] ) +
    _opt_('pull',r'\d+') + _opt_unit_('force',r'\d+','pN') + _opt_time_ + _opt_('series',r'\d+') \
    + _opt_('isBackground',r'background') )

Pattern = re.compile(FILENAME_SYNTAX)

FILENAME_FIELDS = ( 'construct', 'conditions', 'slide', 'mol', 'pull',
  'force', 'min', 'sec', 'series', 'isBackground' )

FileInfo = collections.namedtuple('FileInfo', FILENAME_FIELDS)

# construct, conditions, slide, mol, pull
basePattern = re.compile(r'^(?P<construct>[a-zA-Z0-9]+)_(.*)_s(\d+)(?:m(\d+))?'+_opt_('pull',r'\d+'))

forcePattern = re.compile(r'_(\d+)pN')

# min, sec, series
timePattern = re.compile(r'_(?:(\d+)min)?(?:(\d+)s)?(?:_(\d+))?(?:_|$)')

bgPattern = re.compile(r'_background')

COMMENT_LINE = '#'


def parseFilename(filename):
  #construct,conditions,slide,mol,pull=basePattern.match(filename).groups()
  #slide=int(slide); mol=int(mol); pull=int(pull)
  #force=toNum(forcePattern.search(filename).group(1))
  #min,sec,series=map(toNum, timePattern.search(filename).groups())
  #background = bgPattern.search(filename) is not None

  (construct, conditions, slide, mol, pull, force, min, sec,
    series, isBackground) = Pattern.match(filename).groups()

  slide = toInt(slide)
  mol = toInt(mol)
  pull = toInt(pull) or 1
  force = toNum(force)
  min,sec,series = map(toInt,(min,sec,series))
  series = (series or 1) if min or sec else series
  isBackground = isBackground is not None

  return FileInfo(construct,conditions,slide,mol,pull,force,min,sec,series,isBackground)

LOAD_FUNCTIONS = dict((func.extension,func) for func in 
  filter(lambda x: hasattr(x,'extension'), globals().values()))

