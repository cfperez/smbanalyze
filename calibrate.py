from numpy import loadtxt, array, mean, pi, cos, sin, sqrt, exp, split
from functools import wraps
from scipy.optimize import curve_fit


def load_spectra(filename):
  '''Returns dictionary with all spectra info

Output Format of .spec file:

nmptsr
kxLRO kyLRO
kxVar kyVar
FcX AmpX
dFcX dAmpX
FcY AmpY
dFcY dAmpY
freq powX powY fitX fitY
...
(space delimiters)
  '''
  fit_data_fields = ('kx_lro', 'ky_lro', 'kx_var', 'ky_var', 
      'f_x', 'amp_x', 'df_x', 'damp_x',
      'f_y', 'amp_y', 'df_y', 'damp_y')

  with open(filename) as fh:
      num_pts = int(fh.readline())
      fit_data = []
      for line in fh:
          if line != '' and line.count(' ') > 1:
              break
          fit_data += map(float, line.strip().split())
      return dict(zip(fit_data_fields, fit_data), 
          num_pts=num_pts,
          spectra=loadtxt(fh))

def faxen(a, h):
  """Calculate Faxen's correction to drag
  a = radius in nm
  h = height of bead center in nm
  """
  u = float(a)/h
  return 1/(1-9./16*u + u**3/8. - 45./256*u**4 - u**5/16.)

def position_to_height(position, focal_shift, offset):
  return focal_shift*(position-offset)

def stokes_drag(bead_radius, viscosity):
  "Drag on sphere in SI units. Bead_radius in m, viscosity in cP"
  return 6.*pi*bead_radius*viscosity*1e-3

def stiffness_from_rolloff(rolloff_freq, bead_radius, temperature, height=None):
  '''Returns stiffness in pN/nm calculated from rolloff measurements using faxen's correction at height if given

  bead_radius and height must be in meters
  '''
  faxen_correction = height is None and 1 or faxen(bead_radius, height)
  return 2.*pi*stokes_drag(bead_radius, water_viscosity(temperature))*rolloff_freq*faxen_correction*1e3

def water_viscosity(temp):
  '''Viscosity of water at temp in centiPoise (*1e-3 to get Pa*s)
  '''
  return 0.89e-3 + 2.24e-5*(25.-temp)

def calc_trap_temps(rolloff_strong_on, rolloff_strong_off, 
    power_ratio,
    room_temp=18.33,
    viscosity=water_viscosity):
    '''Calculate temperature increasing due to traps using linear approximation to viscosity around room temp

    power_ratio = WEAK / STRONG

    Returns (exp_T, strong_T, weak_T):
        Experiment Temp
        strong trap delta T
        weak trap delta T

    Solves the equation which relates ratio of rolloff frequencies to inverse ratio of viscosities,
    using a linear approximation to viscosity around room_temp, and assuming the temperature
    scales linearly with the laser power measured before the objective.

'''
    try:
        fbar = 1 - mean(array(rolloff_strong_off)/array(rolloff_strong_on))
    except ValueError:
        raise ValueError('Likely isues with rolloff arrays of unequal size')
    eta_0 = viscosity(room_temp)
    eta_slope = eta_0-viscosity(room_temp+1)
    
    strong_T = eta_0/eta_slope * fbar / (1+power_ratio*fbar)
    weak_T = strong_T * power_ratio
    return strong_T+weak_T+room_temp, strong_T, weak_T

def flatten_once(l):
    return list(level_2 for level_1 in l for level_2 in level_1)

def block(arr, N):
    return array(map(mean, split(arr, N)))

def lp_filter(f3db, order=4):
    f3db = float(f3db)
    return lambda f: 1./(1 + (f/f3db)**(2*order))

def filtered(func, f3db, order=4):
    lp = lp_filter(f3db, order)
    @wraps(func)
    def func_(f, *args, **kwargs):
        return lp(f)*func(f, *args, **kwargs)
    return func_

def P_hydro(f, fc, D, height, viscosity, bead_radius, bead_density):
    f = abs(f)
    l_r = height/bead_radius
    nu = viscosity/1000 # m/s kinematic viscosity
    m_bead = bead_density * 4*pi/3*bead_radius**3
    
    f_v = nu / (pi*bead_radius**2)
    f_m = 6*pi*viscosity*bead_radius / (2*pi * (m_bead + 2*pi*bead_radius**3/3) )
    f_v_ = sqrt(f/f_v)
    Re_gamma = 1 + f_v_ - 3*bead_radius/(16*l_r) + 3/(4*l_r)*exp(-2*l_r*f_v_)*cos(2*l_r*f_v_)
    Im_gamma = -f_v_ + 3/(4*l_r)*exp(-2*l_r*f_v_)*sin(2*l_r*f_v_)
    return Re_gamma/(2*pi**2)*D / ((fc - f*Im_gamma - f**2/f_m)**2 + (f+Re_gamma)**2)

def aliased(func, f_sample, order=3):
    @wraps(func)
    def aliased_func(f, *args):
        return sum(func(f+n*f_sample, *args) for n in xrange(-order,order))
    return aliased_func

def P_hydro_params(param_dict):
    required = set(('height', 'viscosity', 'bead_radius', 
        'bead_density', 'filter_3db', 'f_sample'))
    keys = set(param_dict.keys())
    if set(param_dict.keys()) >= required:
        return param_dict
    else:
        raise ValueError('Missing required parameters: %s' % list(required-keys))

def P_hydro_fit(height, viscosity, bead_radius, bead_density, filter_3db, f_sample, aliasing=False, filter_order=4):
    P_func = filtered(P_hydro, filter_3db, filter_order)
    if aliasing:
        P_alias = aliased(P_func, f_sample)
        return lambda f,fc,D: P_alias(f, fc, D, height, viscosity, bead_radius, bead_density)
    else:
        return lambda f,fc,D: P_func(f, fc, D, height, viscosity, bead_radius, bead_density)
        
def P_lorentz(f, fc, D):
    return D/(2*pi**2)/(fc**2 + f**2)

def fit_power_spectra(frequency, power_, constants, guess=(6000, 1e7)):
    constants = P_hydro_params(constants)
    return curve_fit(P_hydro_fit(**constants), frequency, power_, p0=guess)

def calculate_temperatures_lsq(rolloff_ratio, laser_power_ratio, viscosity=water_viscosity, room_temp=18.33):
    from scipy.optimize import leastsq
    def _temp_min(T):
        low_T = lambda T_: calc_weak_trap_temp(T_, laser_power_ratio, room_temp=room_temp)
        return array(map(
             lambda T_: rolloff_ratio - viscosity(T_)/viscosity(low_T(T_)),
             T
        ))
    out,_ = leastsq(_temp_min, 25)
    high_T = out[0]
    low_T = calc_weak_trap_temp(high_T, laser_power_ratio)
    return high_T, high_T-low_T, low_T-room_temp

def calc_weak_trap_temp(T, power_ratio, room_temp):
    return (T-room_temp)/(1+power_ratio)+room_temp
