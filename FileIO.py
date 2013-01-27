from __future__ import with_statement
import re
import collections
from datetime import datetime
import os
import operator
from inspect import isfunction

import numpy as np

from Types import FretData, PullData
from useful import toNum, dotdict, toInt

IMAGE_FILE = '.img'
CAMERA_FILE = '.cam'
FRET_FILE = '.fret'
PULL_FILE = '.str'

change_extension = lambda x,y: os.path.splitext(x)[0]+y
add_img_ext = lambda x: x+IMAGE_FILE if not x.endswith(IMAGE_FILE) else x
add_cam_ext = lambda x: x+CAMERA_FILE if not x.endswith(CAMERA_FILE) else x
add_fret_ext = lambda x: x+FRET_FILE if not x.endswith(FRET_FILE) else x
add_pull_ext = lambda x: x+PULL_FILE if not x.endswith(PULL_FILE) else x

LOAD_FUNCTIONS = dict((func.extension,func) for func in 
  filter(lambda x: hasattr(x,'extension'), globals().values()))

def load(fname, **kwargs):
  base,extension = os.path.splitext(fname)
  return LOAD_FUNCTIONS[extension](fname,**kwargs)

def loaddat(filename, **kwargs):

  comment_line = kwargs.pop('comments', COMMENT_LINE)
  metaparser = kwargs.pop('metaparser', None)

  colnames = None
  comments = ''

  with open(filename,'rU') as fh:
    position = 0
    for number,line in enumerate(fh):
      if np.any( map(line.startswith, comment_line) ):
        comments += line
      elif line.isspace():
        continue
      elif colnames:
        position = number
        break
      elif any( map(str.isalnum, line.split()) ):
        colnames = line.lower().split()

    fh.seek(0)
    data = np.loadtxt(fh, skiprows=position, **kwargs)

  if metaparser:
    metadata = metaparser(comments) if isfunction(metaparser) else comments
    return metadata, colnames, data

  return colnames, data

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
loadimg.extension=IMAGE_FILE

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


def loadsettings(filename, **kwargs):
  "Return dictionary of key/value pairs from LabView-style configuration file"

  cast = kwargs.get('cast') or str

  settings = dotdict()
  with open(filename) as f:
    header = None
    for line in f:
      line = line.strip()
      if line == '':
        continue

      m = re.match(r'\[(\w+)\]',line)
      if m:
        header = m.group(1).lower()
      else:
        key,value = line.split('=')
        settings[header][key.strip().lower()] = cast(value)

    return settings

def savesettings(filename, file_mode, **settings):
  "save key/value object into LabView-style configuration file"

  with open(filename, file_mode) as f:
    for key in settings:
      f.write( '[%s]\n' % key.capitalize() )

      for key,value in settings[key].iteritems():
        f.write( '%s=%s\n' % (str(key.capitalize()),str(value)) )


def loadfret(fname,**kwargs):
  header,data = loaddat(fname,**kwargs)
  return FretData(*data.T)
loadfret.extension=FRET_FILE

def loadstr(fname,**kwargs):
  header,data = loaddat(fname,comments=('#','/*'),**kwargs)
  if header != ['extension','force','trapdistance']:
    raise IOError, "Stretch file must contain extension, force, and separation"
  return PullData(*data.T)
loadstr.extension=PULL_FILE

##################################################
## Filename Parsing
##################################################
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
