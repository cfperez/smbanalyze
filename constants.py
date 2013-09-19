import logging

parameters = {'T':293.2}

kT = lambda T: 0.0138965*T

beta = 0.13
gamma = 1.16

sumOfBeadRadii = 715
stiffness = (1.525, 0.25)

default_pulling_sampling_time = 0.025
default_fret_exposure_time_ms = 50

default_background_subtract = 94

logHandler = logging.StreamHandler()
logHandler.setFormatter(logging.Formatter('%(asctime)s | %(name)s:%(levelname)s: %(message)s'))
logLevel = logging.WARNING

DEFAULT_FIGURE_EXT = '.png'
