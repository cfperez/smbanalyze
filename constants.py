import logging


temperature = 24.0
parameters = {'T':273.2 + temperature}

kT = lambda T: 0.013806488*T

beta = 0.13
gamma = 1.

sumOfBeadRadii = 710
stiffness = (1.57, 0.275)

default_pulling_sampling_time = 0.01
default_fret_exposure_time_ms = 30

default_background_subtract = 94

logHandler = logging.StreamHandler()
logHandler.setFormatter(logging.Formatter('%(asctime)s | %(name)s:%(levelname)s: %(message)s'))
logLevel = logging.WARNING