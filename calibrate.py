from numpy import loadtxt, array, mean, pi

def load_spectra(filename):
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

def stokes_drag(bead_radius, viscosity):
  "Drag on sphere. Bead_radius in nm, viscosity in cP"
  return 6.*pi*bead_radius*1e-9*viscosity

def stiffness_from_rolloff(rolloff_freq, bead_radius, temperature, height=None):
  faxen_correction = height is None and 1 or faxen(bead_radius, height)
  return 2.*pi*stokes_drag(bead_radius, water_viscosity(temperature))*rolloff_freq*faxen_correction

def water_viscosity(temp):
  return 0.89 + 0.0224*(25.-temp)

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
