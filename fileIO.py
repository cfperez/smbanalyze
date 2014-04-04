import collections
import re
from datetime import datetime
import os
import glob
import ast
from inspect import isfunction

import numpy as np

from useful import toNum, toInt, isInt

IMAGE_FILE = '.img'
CAMERA_FILE = '.cam'
FRET_FILE = '.fret'
PULL_FILE = '.str'
OFC_FILE = '.dat'

REGISTERED_EXT = (IMAGE_FILE,CAMERA_FILE,FRET_FILE,PULL_FILE)

OFC_NUM_OF_COLUMNS = 2

def unique(list_):
  seen = set()
  seen_add = seen.add
  return [x for x in list_ if x not in seen and not seen_add(x)]
uniquify=unique

def sort_files(files):
  return sorted(files, key=lambda f: f.split('_'))

def flist(*globs):
  assert(len(globs)>0)
  return sort_files(glob.glob(makeMatchStrFromArgs(*globs, reg_exp=False)))
filelist = flist

def filter_extensions(files, extensions):
  def condition(s):
    base,ext = os.path.splitext(s)
    return ext in extensions
  return filter(condition, files)
  
def files_matching(globs, with_ext=(), keep_ext=False, unique=False):
  assert(len(globs)>0)
  if not globs:
    globs = ['']
  files = flist(*globs)
  if with_ext:
    files = filter_extensions(files, with_ext)
  if not keep_ext:
    files = [splitext(name)[0] for name in files]
  if unique or not keep_ext:
    files = uniquify(files)
  return files
  
def makeMatchStrFromArgs(*globs, **options):
  reg_exp = options.get('re_match', options.get('reg_exp', True))
  globs = list(globs)
  last = globs[-1]
  if isInt(last):
    globs[-1] = '_'+last
  if reg_exp:
    anychar = '.*'
    endmatch = r'(_|$)'
  else:
    anychar = '*'
    endmatch = '*'
  return r'{any}{pattern}{end}'.format(pattern=anychar.join(globs), any=anychar, end=endmatch)


class fileIOError(Exception):
  IGNORE = 2
  WARNING = 1
  ERROR = 0

  def __init__(self, message, level=ERROR):
    self.level = level
    self.strerror = message

  @property
  def isError(self):
    return self.level == fileIOError.ERROR

def toSettings(comments):
  # Careful! Using literal_eval to cast so that lines can contain
  # lists e.g. stiffness = [1,1.5]
  return parseSettingsFromLines(comments.splitlines(), parse_config_string)

def string_is_list(str_):
  return str_.startswith('(') or str_.startswith('[')
  
def parse_config_string(value):
  if not isinstance(value, str):
    return value
  if string_is_list(value):
    value_list = value.strip('()[]').split(',')
    return map(toNum, value_list)
  elif value:
    return ast.literal_eval(value)
  else:
    return value

def _to_datetime(date_list):
  if date_list[0] < 2000:
    date_list[0] += 2000
  return datetime(*date_list)

def load(fname, comments=toSettings, header=False, **kwargs):
  if not os.path.exists(fname):
    raise IOError('No file named "{}" found'.format(fname))
  base,extension = os.path.splitext(fname)
  fromLoad =  LOAD_FUNCTIONS[extension](fname, **kwargs)
  if isinstance(fromLoad, tuple):
    loadComments, loadHeader, loadData = fromLoad
    try:
      loadComments = comments(loadComments) if isfunction(comments) \
                      else loadComments
    except ValueError:
      raise fileIOError(level=fileIOError.WARNING, 
          message='Comments from file {0} could not be processed'.format(fname))
    output = [loadData]
    if header:
      output.insert(0,loadHeader)
    if comments:
      try:
        if 'date' in loadComments:
          loadComments['date'] = _to_datetime(loadComments['date'])
      except e:
        print "WARNING parsing comments: %s" % e
      output.insert(0,loadComments)
    return tuple(output) if len(output)>1 else output[0]
  else:
    return fromLoad

def loadstr(fname, **loadOptions):
  loadOptions.setdefault('comments', ('/*', '#'))
  filecomments,fileheader,data = loaddat(fname, **loadOptions)
  if fileheader != ['extension','force','trapdistance']:
    raise fileIOError("Stretch file must contain extension, force, and separation")
  return filecomments, fileheader, data
loadstr.extension=PULL_FILE

def loaddat(filename, **kwargs):

  comment_line = kwargs.pop('comments', COMMENT_LINE)

  colnames = None
  comments = ''
  with open(filename,'rU') as fh:
    position = 0
    # loop stops after grabbing column names
    for number,line in enumerate(fh):
      if np.any( map(line.startswith, comment_line) ):
        if not line.isspace():
          comments += line.strip(''.join(comment_line))
      else:
        line = line.replace('\\t',' ')
        if any( map(str.isalnum, line.split()) ):
          colnames = line.lower().split()
          position = number+1
          break

    fh.seek(0)
    data = np.loadtxt(fh, skiprows=position, **kwargs)

  return comments, colnames, data

def loadfret(filename, **kwargs):
  comments,header,data = loaddat(filename,**kwargs)
  return comments, header, data
loadfret.extension=FRET_FILE

def savefret(filename, time, donor, acceptor, fret, metadata=None, comments=''):
  if comments:
    comments += '\n'
  if metadata:
    comments += fromSettings(metadata)
  savedat(filename, (time,donor,acceptor,fret),
    header='time donor acceptor FRET', 
    fmt=('%.3f','%u','%u','%.5f'),
    comments=comments)

def loadsettings(filename, **kwargs):
  "Return dictionary of key/value pairs from LabView-style configuration file"
  cast = kwargs.get('cast') or str
  with open(filename, 'rU') as f:
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
    return header, None, None
  elif '=' in line:
    key,value = line.split('=')
    return None, key.strip().lower().replace(' ','_'),value.strip()
  else:
    return None, line, None

def strip_blank_and_comments(iter_str):
  for line in strip_blank(iter_str):
    if not line.startswith('#'): yield line

def strip_blank(iter_str):
  for line in iter_str:
    line=line.strip()
    if line != '': yield line

def savedat(filename, data, header='', comments='', fmt='%.9e', delimiter='\t'):

  newline = '\n' if comments else ''
  header = ''.join(map(lambda line: COMMENT_LINE[0]+' '+line, comments.splitlines(True))) \
      + newline + header.replace(' ',delimiter)
  if isinstance(data, tuple):
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
    raise fileIOError("Image file %s is corrupted, expected frame: %d, height: %d, width: %d" % 
        (filename, img_size[0], img_size[1], img_size[2]))
  return img

def loadofc(filename, datatype='>f4', **kwargs):
  data = np.fromfile(filename, datatype)
  num_columns = OFC_NUM_OF_COLUMNS
  try:
    data_reshaped = data.reshape( (-1,num_columns) )
  except ValueError:
    raise fileIOError(
      "OFC file is wrong format or is corrupted; expect {} columns".format(num_columns))
  return data_reshaped

def loadcam(filename):
  "Return dictionary of camera settings as contained in .cam file"
  # Can be replaced by loadsettings and a handle for the DateTime key value
  # once the old-style .cam setting is supplanted

  settings = {}
  with open(filename, 'rU') as f:
    for line in f.readlines():
      key,value = line.strip().split('\t')
      m = re.match('^\d+/', key)
      if key.lower() != 'datetime' and not m:
        value = int(value)
      else:
        if m: # Handle old style cam files which don't have a DateTime key
          value = key + '\t' + value
          key = 'datetime'
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
        if setting == 'datetime': continue
        output += SETTING_FMT % (setting,value)
    except AttributeError:
      # UPDATE ME!
      if header == 'datetime': continue
      output += SETTING_FMT % (header, subHead)
  return ''.join(output)

##################################################
## Filename Parsing
##################################################
def splitext(fname):
  if not fname:
    return None, None
  fname = os.path.basename(fname)
  basename,ext=os.path.splitext(fname)
  if ext not in REGISTERED_EXT:
    basename,ext=fname,''
  return basename,ext

change_extension = lambda x,y: os.path.splitext(x)[0]+y
add_img_ext = lambda x: x+IMAGE_FILE if not x.endswith(IMAGE_FILE) else x
add_cam_ext = lambda x: x+CAMERA_FILE if not x.endswith(CAMERA_FILE) else x
add_fret_ext = lambda x: x+FRET_FILE if not x.endswith(FRET_FILE) else x
add_pull_ext = lambda x: x+PULL_FILE if not x.endswith(PULL_FILE) else x

def filesFromName(filename, *extensions):
  '''Return basename + filenames of given extensions (including original) using basename from filename.
  Default is to return PULL_FILE and FRET_FILE extensions.'''
  assert isinstance(filename, str)
  basename,ext = splitext(filename)
  if basename == '':
    raise ValueError('Filename "{}" must contain a parseable basename'.format(filename))

  if len(extensions) == 0:
    extensions = (PULL_FILE, FRET_FILE)
  if ext and ext not in extensions:
    extensions = (ext,) + extensions

  return [basename] + [basename+ext for ext in extensions]


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
    _opt_('pull',r'\d+') + _opt_unit_('force',r'\d+(?:.\d+)?','pN') + _opt_time_ + _opt_('series',r'\d+') \
    + _opt_('isBackground',r'background') )

Pattern = re.compile(FILENAME_SYNTAX)

FILENAME_FIELDS = ( 'construct', 'conditions', 'slide', 'mol', 'pull',
  'force', 'min', 'sec', 'series', 'isBackground' )

FileInfo = collections.namedtuple('FileInfo', FILENAME_FIELDS)

# construct, conditions, slide, mol, pull
basePattern = re.compile(r'^(?P<construct>[a-zA-Z0-9]+)_(.*)_s(\d+)(?:m(\d+))?'+_opt_('pull',r'\d+'))

forcePattern = re.compile(r'_(\d+(?:\.\d+)?)pN')

# min, sec, series
timePattern = re.compile(r'_(?:(\d+)min)?(?:(\d+)s)?(?:_(\d+))?(?:_|$)')

bgPattern = re.compile(r'_background')

COMMENT_LINE = ('#','/*')

MOL_FILENAME_INFO = ('construct', 'conditions', 'slide', 'mol')
MOL_FILENAME_NUM_FIELDS = 3

def split_on_token(string, token='_', default_not_found=None, maxsplit=None):
  split_args = (token, maxsplit) if maxsplit else (token,)
  split = string.split(*split_args) if token in string else (string, default_not_found)
  if not maxsplit:
    return split
  return split + [default_not_found] * (1+maxsplit-len(split))

def make_info_dict(values):
  assert len(MOL_FILENAME_INFO) == len(values)
  return dict(zip(MOL_FILENAME_INFO, values))

mol_info_match = re.compile(r'(?P<info>.*s\d+m\d+)_?(?P<details>.*)').match

def parse_mol_info2(filename):
  fields = split_on_token(filename, maxsplit=MOL_FILENAME_NUM_FIELDS)
  construct, condition, slidemol, the_rest = fields
  slide, mol = map(int, (slidemol[1], slidemol[3]))
  return make_info_dict((construct, condition, slide, mol)), the_rest

def parse_slide_mol_str_to_ints(slide_mol):
  return map(int, slide_mol.lower()[1:].split('m'))

def split_outside_in(string_, delimiter='_'):
  '''a_b_d_e_f => (a, [b,d,e], f)'''
  string_split = string_.split(delimiter)
  return string_split[0], string_split[1:-1], string_split[-1]

def parse_mol_info(filename):
  mol_info_str, the_rest = mol_info_match(filename).groups()
  construct, conditions, slidemol = split_outside_in(mol_info_str)
  slide, mol = map(int, slidemol.lower()[1:].split('m')) #parse_slide_mol_str_to_ints(slidemol)
  return make_info_dict((construct, conditions, slide, mol)), the_rest

def parseFilename(filename):
  if not filename:
    return None
  try:
    (construct, conditions, slide, mol, pull, force, min, sec,
      series, isBackground) = Pattern.match(filename).groups()
  except AttributeError:
    return None

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

