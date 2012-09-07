import Image, FileIO
import matplotlib.pyplot as plt
from FileIO import savedat

beta = 0.13
gamma = 1.16

def calc(stack, beta=beta, gamma=gamma):
    """Calculates FRET of a pull from an Image.Stack

calcFRET( Image.Stack, beta = Image.beta, gamma = Image.gamma)

defaults: taken from previous measurements

RETURNS array of calculated FRET for each frame
"""

    donor = stack.donor - min(stack.donor)
    acceptor = stack.acceptor - donor*beta
    acceptor = acceptor - min(acceptor)

    return acceptor/(acceptor+gamma*donor)

def tofile(stack, filename, **kwargs):
    "saveFRETdata( fret, ImageStack, filename): saves donor,acceptor, FRET to 3 column text file"

    fret = pullfromStack(stack, **kwargs)
    savedat(filename, (stack.donor,stack.acceptor,fret), header='donor acceptor FRET', fmt=('%u','%u','%.5f'))
    return fret
