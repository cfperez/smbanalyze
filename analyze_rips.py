# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:43:38 2013

@author: Christian
"""
import os
import os.path as path
import sys
from contextlib import contextmanager
from smbanalyze import experiment, browser
from numpy import r_, mean, median, std, array
import matplotlib.pyplot as plt

WORKING_DIR = '/Volumes/users2/Force-Fret Project/Analysis/!Rip Analysis/'

# curve adjustment parameters
EXTENSION_OFFSET = 1050

# For outliers
MIN_FORCE_CUTOFF = 8.5
MAX_FORCE_CUTOFF = 14

# For FEC fits
HSTART = 900
HEND = 1010
RSTART = 1035
REND_FORCE = 15

def calc_rip_size(pulls):
    p = pulls.collapse()
    p.fitHandles(x=(HSTART, HEND))
    return p.fitRip(RSTART, REND_FORCE)['Lc1']

def calc_rip_forces(pulls):
    return pulls.findRip(HEND)[:,1].tolist()

def process_experiments(pulls):
    'Return rip size and rip forces from list of ExperimentBrowser'
    assert isinstance(pulls, experiment.List)
    print "Calculating {}".format(pulls[0].filename)
    pulls.adjustOffset(to_x=EXTENSION_OFFSET)
    return calc_rip_size(pulls), calc_rip_forces(pulls)

# Close the figure used by plotall() if the variable "output" exists
if "output" in dir() and hasattr(output, 'figure'):
    output.figure.close()
    
output = experiment.List()
rip_size = []
rip_forces = []

filename = path.join(WORKING_DIR, sys.argv[1])
with browser.ExperimentBrowser(filename) as exp_browser:
    for mol in exp_browser:
        mol = mol.has_value(trap_f_atleast=16)
        mol.adjustOffset(to_x=EXTENSION_OFFSET)
        output += mol
        rip_size += [calc_rip_size(mol)]
        rip_forces += calc_rip_forces(mol)

with ChangeDirectory(WORKING_DIR):
    plt.clf()
    output.plotall(legend=None)
    figure_name = path.splitext(filename)[0]
    output.savefig(figure_name+'_pulls')

    plt.figure()
    plt.clf()
    plt.hist(rip_forces, bins=14, range=(MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
    plt.title('Rip force distribution')
    plt.show()
    plt.xlim((MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
    plt.savefig(figure_name+'_rip_force_histogram')

    plt.figure()
    plt.hist(rip_size)
    plt.show()

    print "\nCalculating stats from {} traces".format(len(rip_forces))
    print "Mean rip stats:\n\tforce = {:.2f} ({:.2f})\trip size = {:.2f} ({:.2f})".format(
        mean(rip_forces), std(rip_forces), 
        mean(rip_size), std(rip_size))
    print "Median rip stats:\n\tforce = {:.2f}\tsize = {:.2f}".format(
        median(rip_forces), median(rip_size))