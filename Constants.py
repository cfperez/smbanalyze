import logging

parameters = {'T':293.2}

kT = lambda T: 0.0138965*T

beta = 0.13
gamma = 1.16

default_background_subtract = 100

logHandler = logging.StreamHandler()
logHandler.setFormatter(logging.Formatter('%(asctime)s | %(name)s:%(levelname)s: %(message)s'))
logLevel = logging.DEBUG
