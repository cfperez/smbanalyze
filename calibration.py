from __future__ import with_statement
import numpy.core.records as rec
from numpy import array, average, sort, hstack, reshape, linspace, loadtxt, savetxt, iterable, float, pi, abs
import scipy.interpolate as interpolate
from pylab import plot, figure, subplot, legend, xlabel, ylabel
from matplotlib import lines
from scipy import polyfit, poly1d
import os
from curvefit import Fit

__all__ = ['CalibrateTrap', 'average_waves']

STD_TEMPLATE = {'lro':'KvsAmp%d.dat', 'var':'VarKvsAmp%d.dat', \
	'stokes': 'StokesDrag_bead%d_gain%d.dat'}

FOCAL_SHIFT = 0.82

def faxen(a, h):
  """ Calculate Faxen's correction to drag
  a = radius in nm
  h = height of bead center in nm
  """
  u = float(a)/h
  return 1/(1-9./16*u + u**3/8. - 45./256*u**4 - u**5/16.)

def freq_vs_height(stage, F_actual=5000., surface=10000, bead_radius=300):
  "FITFUNC: stage (nm), F_actual (Hz), surface (nm), bead_radius (nm)"
  return F_actual/faxen(bead_radius, FOCAL_SHIFT*(surface-stage)+bead_radius)

def position_to_height(position, focal_shift, offset):
  return focal_shift*(position-offset)

def freq_vs_position(position, F_at_infinity=5000, focal_shift=FOCAL_SHIFT, height_offset=150, bead_radius=300):
  "FITFUNC: position (nm), F_actual (Hz)"
  return F_at_infinity/faxen(bead_radius, position_to_height(position, focal_shift, height_offset))

def fit_vs_height(height, Y, fitfunc, **fit_params):
  fit = Fit(height, Y, fitfunc,
      fixed='bead_radius',
      **fit_params)
  return fit

def fit_freq_vs_position(position, frequency, **fit_params):
  default_fixed = ('focal_shift', 'bead_radius')
  fixed = fit_params.pop('fixed', default_fixed)
  return Fit(position, frequency, freq_vs_position,
    fixed=fixed,
    **fit_params)

def fit_freq_vs_height(stage, freq, bead_radius, **fit_params):
  "Return rolloff frequency from Faxen vs height calibration of single bead"
  fit = Fit(stage, freq, freq_vs_height,
    fixed='bead_radius', bead_radius=bead_radius,
    **fit_params)
  return fit

def drag(bead_radius, viscosity):
  "Drag on sphere. Bead_radius in nm, viscosity in cP"
  return 6.*pi*bead_radius*1e-9*viscosity

def water_viscosity(temp):
  return 0.89 + 0.0224*(25.-temp)

def stiffness_from_rolloff(rolloff_freq, bead_radius, temperature, height=None):
  faxen_correction = height is None and 1 or faxen(bead_radius, height)
  return 2.*pi*drag(bead_radius, water_viscosity(temperature))*rolloff_freq*faxen_correction

class TrapCalibration(object):
  def __init__(self, bead_diameter, temp, focal_shift=FOCAL_SHIFT):
    ''' bead_diameter in nm'''
    self.bead_radius = bead_diameter / 2.
    self.temp = temp
    self.viscosity = water_viscosity(temp)
    self.focal_shift = focal_shift
    self.stiffness = []

  def fromFreqVsHeight(self, height, freq, **fit_params):
    self.fit = fit_freq_vs_height(height, freq, self.bead_radius, **fit_params)
    self.rolloff = self.fit['F_actual']
    self.stiffness.append(
        stiffness_from_rolloff(self.rolloff, self.bead_radius, self.viscosity)
    )
    return self.stiffness[-1]

  def from_freq_vs_position(self, position, freq, **starting_params):
    self.fit = fit_freq_vs_position(height, freq, 
        bead_radius=self.bead_radius, 
        focal_shift=self.focal_shift,
        **starting_params)
    self.rolloff = self.fit['F_actual']
    self.stiffness.append(
        stiffness_from_rolloff(self.rolloff, self.bead_radius, self.viscosity)
    )
    return self.stiffness[-1]

  def average_stiffness(self):
    return average(self.stiffness)

class CalibrateTrap:
    """
    Description:

	Loads calibration data from LRO, Variance, and Stokes measurements
	and calculates the average stiffness versus gain, and fits this stiffness
	to the measured laser power.

	To load data from multiple files, the class uses a template system that
	mirrors the standard naming scheme of the VIs:

	    KvsAmp%d.dat
	    VarKvsAmp%d.dat
	    StokesDrag_bead%d_gain%d.dat

	All of these can be changed by passing in named parameter to the class
	constructor containing keys lro, var, or stokes and the desired template
	string:
	    
	    >>> trap = CalibrateTrap( 0.73, stokes='bead%d_gain%d.dat' )

    Methods:

    __init__(bead_size, path='', **file_prefix)

	bead_size: in microns

    lro(file_list, output='LRO_avg.dat', template=None)
    var(file_list, output='var_avg.dat', template=None)
    stokes(bead, gain, template=None)

	Loads data from various dat files. You can override the default template here.

    load_all(**kwargs)

	Loads data from all files in one go, using keyword arguments to supply the
	file_list, bead and gain parameters to the three functions. Requires use of
	standard templates or templates updated at object instantiation.

	    load_all(lro=file_list, var=file_list, stokes=(bead, gain))
	    load_all(lrovar=file_list, stokes=(bead,gain))

    calculate(faxen)

	Apply faxen's correction to measurements and find average stiffness vs gain.

	faxen: faxen's correction for the given trap

	Returns k_avg (2D array for Kx and Ky in same order as self.gain)

    fitKvsPower(filename="PowerVGain.dat", zeroed="_0")

	Load laser power vs gain file, background subtract by using filename with 
	zeroed appended to filename, and return fit parameters for kx and ky
	versus power.

    save(filename='calibration.dat', delimiter='\t')

	Save all stiffness, gain, and power measurements in an Igor compatible file

    Attributes:

	path: current path

	templates: string templates for file loading

	bead_size: self.explanatory

	LRO_avg: averaged measurements using LRO method

	var_avg: averaged measurements using Var method

	stokes_avg: averaged stiffnesses using Stokes method

	stokes_inter: function to allow interpolation of Stokes stiffness over measured
			range

	AFTER running calculate()
	-------------------------

	gain: gains of each AOD calibrated over
	k_avg: stiffnesses for X and Y corresponding to gains
    """
    def __init__(self, bead_size, path='', **file_prefix):
        #self.templates = file_prefix or STD_TEMPLATE
	self.templates = STD_TEMPLATE
	self.templates.update(file_prefix)
	self.path = path
	self.gain = []
	self.bead_size = bead_size

    def lro(self, file_list, output='LRO_avg.dat', template=None):
	"""
	lro(file_list, output='LRO_avg.dat', template=None)

	Loads data from various lro dat files. You can override the default template.

	file_list: either a list of numbers to apply to the template, or a single number that
	the function will automatically count up to (i.e. 6 will load files 1 through 6)

	Saves averaged data to filename given by output
	"""

	template = template or self.templates['lro']
	if type(file_list) != list:
	    file_list = fullrange(file_list)[1:]

	self.__lro_file_list = file_list

	files = [ self.path + template % x for x in file_list ]
	print "Loading files: " + str(files)
	self.LRO_avg, self.LRO_raw = average_waves(files, output=self.path + output, usecols=(0,3,4))
	self.LROvar_avg, self.LROvar_raw = average_waves(files, output=self.path+output.replace('LRO','LROvar'), usecols=(0,1,2))

    def var(self, file_list, output='var_avg.dat', template=None):
	"""
	var(file_list, output='var_avg.dat', template=None)

	Loads data from var dat files. You can override the default template.

	file_list: either a list of numbers to apply to the template, or a single number that
	the function will automatically count up to (i.e. 6 will load files 1 through 6)

	Saves averaged data to filename given by output
	"""
	template = template or self.templates['var']
	if type(file_list) != list:
	    file_list = fullrange(file_list)[1:]

	self.__var_file_list = file_list

	files = [ self.path + template % x for x in file_list ]
	print "Loading files: " + str(files)
	self.var_avg, self.var_raw = average_waves(files, output=self.path + output, usecols=[0,1,2])

    def stokes(self, bead, gain, output='stokes_avg.dat', template=None, path='Stokes'):
	"""
	stokes(bead, gain, output='stokes_avg.dat', template=None, path='Stokes')

	Loads data from stokes dat files. You can override the default template.

	bead & gain: similar to file_list in usage, but specificly fills out the first and
	second template positions for stokes drag templates

	Saves averaged data to filename given by output
	"""
	template = template or self.templates['stokes']

	#if type(bead) != list:
	if not iterable(bead):
	    bead = fullrange(bead)[1:]
	#if type(gain) != list:
	if not iterable(gain):
	    gain = fullrange(gain)[1:]

	stokes_k = []
	stokes_gain = []
	for g in gain:
	    files = [ os.path.join(self.path, path, template) % (b,g) for b in bead ]
	    print "Loading files: " + str(files), "\n"

	    #---------------------------------
	    # Load data from each Stokes file
	    # and calculate stiffness. Store
	    # values in object variables
	    #---------------------------------
	    k = []
	    for file in files:
		try:
		    with open(file, 'r') as f:
			header = [ f.readline() for dummy in range(9) ]
			vel, x, y = loadtxt(f, unpack=True)
		except IOError, (errno, errstr):
		    if errno == 2:
			print "File %s not found! Skipping..." % file
		    else:
			raise
		else:
		    pfit = polyfit(vel*1000, y, 1)
		    k.append(self.bead_size * 9.42e-6 / pfit[0])
		    print file, "=", str(poly1d(pfit, variable='v')).strip(), "stiffness = %.3f" % k[-1]

	    if k == []:
		continue

	    k = -average(k)
	    stokes_k.append(k)
	    stokes_gain += [ int(header[4].split()[2]) ]
	    print "\nAverage k for gain %.1f: %f" % (stokes_gain[-1], k), "\n--------\n"

	# Outside for loop

	# create k vs gain for stokes measurements and sort by gain
	# (interpolation function fails unless gain is strictly increasing)
	self.stokes_avg = sort(rec.fromarrays([stokes_gain, stokes_k], names='gain,k'), order='gain')

	self.stokes_inter = interpolate.interp1d(self.stokes_avg['gain'], self.stokes_avg['k'])
	newX = linspace( self.stokes_avg['gain'][0], self.stokes_avg['gain'][-1], num=15 )
	if hasattr(self, 'LRO_avg'):
	    newX = self.LRO_avg[:,0]
	elif hasattr(self, 'var_avg'):
	    newX = self.var_avg[:,0]

	savetxt(self.path + output, self.stokes_inter(newX) )

    #============================================================================================
    def fitKvsPower(self, filename="PowerVGain.dat", zeroed="_0"):
	
	if not hasattr(self, 'k_avg'):
	    raise RuntimeError, "Must load all data and run calculate() before fitting"

	gain, self.power = loadtxt( filename, skiprows=1, unpack=True )

	if zeroed != None:
	    gain_0, power_0 = loadtxt( filename.replace('.', zeroed + '.', 1), skiprows=1, unpack=True )
	    self.power = self.power - power_0
	    savetxt(filename.replace('.', '_subtract.'), self.power, fmt='%.5f')

	if self.gain.tolist() != gain.tolist():
	    raise ValueError, "Laser power was not measured over the same gains as calibration measurements:\n\
		Laser gain: " + str(gain) + "\nMeasurement gains: " + str(self.gain)

	self.kfitX = polyfit(self.power, self.k_avg[:,0], 1)
	self.kfitY = polyfit(self.power, self.k_avg[:,1], 1)
	print "Fit parameters for X axis:"
	print poly1d(self.kfitX, variable='power')
	print "\nFit parameters for Y axis:"
	print poly1d(self.kfitY, variable='power')


    #============================================================================================
    def calculate(self, faxen):
	"""Calculate the average stiffness using all input files.

	Usage: calculate(faxen)
	
	Description:

	After verifying LRO, variance, and stokes measurements have been loaded,
	this function applies faxen's correction and calculates the average stiffness

	Parameters:

	    faxen

		Faxen's correction for the given trap, depends on bead size
	"""
	try:
	    self.LRO_avg
	    self.LROvar_avg
	    self.var_avg
	    #self.stokes_avg
	except NameError:
	    raise RuntimeError, "Must run lro and var (stokes optional) methods before calculating stiffness"

	# check for proper gains on all files
	if self.LRO_avg[:,0].tolist() != self.var_avg[:,0].tolist():
	    raise ValueError, "File from {0} and {1} do not record the same gains".format(\
		self.templates['lro'], self.templates['var'] ) + \
		"\n\nLRO gains {0}\nVar gains {1}".format(\
		self.LRO_avg[:,0], self.var_avg[:,0])

	self.gain = self.LRO_avg[:,0]

	bead_size = self.bead_size

	# Keep only variance stiffness data
	self.varK_avg = self.var_avg[:,1:3]

	# average variance of two variance techniques
	#self.variance_k_avg = average([self.LRO_avg[:,1:3], self.var_avg[:,:2]], axis=0)
	self.variance_k_avg = ( self.LROvar_avg[:,1:] + self.varK_avg ) / 2

	# apply faxen's correction to roll off data
	self.LROK_avg = self.LRO_avg[:,1:] / faxen

	# combine roll off and variance. Average taken after stokes is accounted for
	self.k_avg = self.LROK_avg + self.variance_k_avg # includes x and y

	if hasattr(self, 'stokes_avg'):
	    try:
		# interpolate gains for stokes using same gains as other methods
		stokes_k_interpolated = self.stokes_inter( self.gain )
	    except ValueError:
		raise ValueError, "Stokes measurements were not done over same gain interval as other measurements!"

	    stokes_k_interpolated /= faxen
	    self.stokes_avg = stokes_k_interpolated

	    self.k_avg[:,1] += stokes_k_interpolated # only includes y axis

	    # list of stiffness values for X and Y
	    # using average of 2 methods for X and 3 for Y
	    self.k_avg /= [2,3]
	else:
	    # list of stiffness values for X and Y
	    # using average of 2 methods for X and Y (no stokes)
	    self.k_avg /= 2
	    

	return hstack((self.gain.reshape(len(self.k_avg),1),self.k_avg))


    def load_all(self, **arg_dict):
	try:
	    self.stokes( *arg_dict['stokes'] )
	    self.lro( arg_dict['lro'] )
	    self.var( arg_dict['var'] )
	except KeyError, (errno, errstr):
	    if 'lrovar' in arg_dict and errstr == 'lro':
		self.lro( arg_dict['lrovar'] )
		self.var( arg_dict['lrovar'] )
	    else:
		raise KeyError, "Missing argument list for %s" % errstr

    def save(self, filename='calibration.dat', delimiter='\t', raw=True):

	# save calculated values in comments
	comments=\
	    """# ky vs power: %s # kx vs power: %s\n""" % \
		(str(poly1d(self.kfitX, variable='power')).strip(), 
		str(poly1d(self.kfitY, variable='power')).strip() )

	# patterns for col names of raw data
	lro_raw_headings = ['kxVar%d', 'kyVar%d', 'kxLRO%d', 'kyLRO%d']
	var_raw_headings = ['kx%d', 'ky%d']

	# col names for raw data
	if raw:
	    apply_to = lambda mylist, value: [ x % value for x in mylist ]
	    lro_head = []
	    var_head = []
	    for lronum, varnum in zip(self.__lro_file_list, self.__var_file_list):
		lro_head += apply_to(lro_raw_headings, lronum)
		var_head += apply_to(var_raw_headings, varnum)
	

	# col names for average values
	header = [ 'gain', 'power', 'kx_avg', 'ky_avg', 'kxLRO_avg', 'kyLRO_avg', 'kxvar_avg', 'kyvar_avg', 'ky_stokes' ]
	if raw:
	    header += lro_head + var_head

	# Tuple of data to be saved in the same order as header listing
	data = (self.gain, self.power, self.k_avg, self.LRO_avg, self.variance_k_avg, self.stokes_avg)

	if raw:
	    # Except first column (gain), append raw data loaded from each file
	    data += (self.LRO_raw[:,1:], self.var_raw[:,1:])

	# data variable must be reshaped so it can be stacked column ways when saved
	def my_shape(x):
	    if x.ndim > 1:
		return x
	    else:
		return reshape(x, (-1, 1))

	save_data = hstack(map(my_shape, data))

	with open(filename, 'w') as fh:
	    fh.write(comments)
	    fh.write( delimiter.join(header) + '\n')
	    savetxt( fh, save_data, fmt='%.6e', delimiter=delimiter)

	print "Data saved to file %s" % filename
	
    def show(self, **to_plot):

	marker = ('+', 's', 'o', '', 'x')
	colors = ('c', 'g', 'k', 'k', 'r')
	names = ('lro', 'var', 'avg', 'fit', 'stokes')

	line_style = map(lambda x,y: x+y, colors, marker)
	styles = dict( zip(names, line_style) )

	for key in to_plot:
	    if not styles.has_key(key):
		raise ValueError, "show() can only be asked to display data for %s" % ', '.join(names)

	styles.update(to_plot)
	#styles looks like {'lro': 'c+', 'var': 'gs', 'stokes': 'rx', 'avg': 'ko', 'fit': 'k' }

	figure()

	# Plot x axis stiffness vs power
	subplot(211)

	lineX = plot(self.power, self.LRO_avg[:,0], styles['lro'], self.power, self.variance_k_avg[:,0], styles['var'], \
	    self.power, self.k_avg[:,0], styles['avg'], self.power, poly1d(self.kfitX)(self.power), styles['fit'])

	map( lines.Line2D.set_label, lineX, names[:-1] )
	map( lines.Line2D.set_markersize, lineX, [10] * len(names[:-1]) )

	legend(loc='upper left')
	xlabel('Power (arb units)')
	ylabel('X stiffness (pN/nm)')

	# Plot y axis stiffness vs power
	subplot(212)
	lineY = plot(self.power, self.LRO_avg[:,1], styles['lro'], self.power, self.variance_k_avg[:,1], styles['var'], \
	    self.power, self.k_avg[:,1], styles['avg'], self.power, poly1d(self.kfitY)(self.power), styles['fit'], \
	    self.power, self.stokes_avg, styles['stokes'] )

	# apply line labels
	map( lines.Line2D.set_label, lineY, names )
	# up marker size to 2
	map( lines.Line2D.set_markersize, lineY, [10] * len(names) )

	legend(loc='upper left')
	xlabel('Power (arb units)')
	ylabel('Y stiffness (pN/nm)')

	return lineY

#---------------------------------------------------------------------------------
def fullrange(*args):
    if len(args) == 1:
	stop = [args[0]]
    else:
	stop = [args[1]]
	if len(args)==3 and (stop[0]-args[0]) % args[2] != 0:
	    stop = []

    return range(*args) + stop

def average_waves(files, output=None, skip=1, usecols=None):
    "Average .dat files together and write to file named output or files.prefix+Avg+files.base"
    
    AvgData = array([])
    AllData = array([])
    avg_over = len(files)
    for file in files:
        
        try:
	    data = loadtxt(file, skiprows=skip, usecols=usecols)
            AvgData = AvgData + data
	except IOError:
	    print "File not found: {0}. Skipping...".format(file)
	    avg_over -= 1
	    continue
        except ValueError:  #   This happens when the two files are not the same size; data is skipped in averaging
            if AvgData.size == 0:
                AvgData = data
            else:
                print "Data file is wrong size: skipping %s" % file
                avg_over -= -1
		continue
	except TypeError:
	    # when there's text in the file
	    avg_over -= 1

	if AllData.size != 0:
	    AllData = hstack((AllData,data[:,1:]))
	else:
	    AllData = data
            
    AvgData = AvgData / avg_over
    
    if output:
        print "Writing averaged waves to file %s" % output
        savetxt(output, AvgData, fmt='%.6E')
    
    return AvgData, AllData
