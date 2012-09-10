from __future__ import with_statement
from useful import toNum, dotdict
import numpy as np
from datetime import datetime
import re

_opt_ = lambda p: r'(?:_(%s))?' % p
__build_pattern__ = lambda *x: r''.join(x)+'$'

FILENAME_SYNTAX = __build_pattern__( \
    r'^([a-zA-Z0-9]+)_.*s(\d+)m(\d+)' + _opt_(r'\d+') + _opt_(r'\d+min') + \
    _opt_(r'\d+') + _opt_(r'background') )


def loadat(filename, comments='#', **kwargs):

    colnames = None
    with open(filename, 'r') as fh:
	for line in fh:
	    if line.startswith(comments) or line.isspace():
		continue

	    if not colnames and any( map(str.isalpha, line.split()) ):
		colnames = line.lower().split()
	    else:
		break

	data = load(fh, **kwargs)
    # end with

    #data.dtype = dtype( zip(colnames, ['f8']*len(colnames)) )

    return data

def savedat(filename, data, header='', comments='#', fmt='%.18e', delimiter='\t'):

    header = ''.join(map(lambda line: comments+line.replace(' ',delimiter), header.splitlines(True)))

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
	raise IOError, "Image file %s is corrupted, expected frame: %d, height: %d, width: %d" % \
	    (filename, img_size[0], img_size[1], img_size[2])
    
    return img

def loadcam(filename):
    "Return dictionary of camera settings as contained in .cam file"

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


def loadsettings(filename, **kwargs):
    "Load key/value objects from LabView-style configuration file"

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
		if not settings.has_key(header):
		    settings[header] = dotdict()
	    else:
		key,value = line.split('=')
		settings[header][key.lower()] = cast(value)

    return settings

def savesettings(filename, **settings):
    "save key/value object into LabView-style configuration file"

    file_mode = settings.get('file_mode','w')

    with open(filename, file_mode) as f:
	for key in settings:
	    f.write( '[%s]\n' % key.capitalize() )

	    for key,value in settings[key].iteritems():
		f.write( '%s=%s\n' % (str(key.capitalize()),str(value)) )
