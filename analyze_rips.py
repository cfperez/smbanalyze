# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:43:38 2013

@author: Christian
"""
import os
import os.path as path
import sys
from contextlib import contextmanager
from smbanalyze import experiment
from numpy import r_, mean, median, std
import matplotlib.pyplot as plt

# curve adjustment parameters
EXTENSION_OFFSET = 1050

# For outliers
MIN_FORCE_CUTOFF = 8.5
MAX_FORCE_CUTOFF = 14

# For FEC fits
HSTART = 900
HEND = 1015
RSTART = 1035
REND_FORCE = 15

# Change the drive letter for your machine!
DATA_DIR = '/Volumes/users2/Force-FRET Project/Data'
ANALYSIS_DIR = '/Volumes/users2/Force-FRET Project/Analysis/!Rip Analysis'

pulls_filename = path.join(ANALYSIS_DIR, sys.argv[1])
with open(pulls_filename) as fh:
    good_mol = [line.strip().split() for line in fh.readlines()]

@contextmanager
def ChangeDirectory(directory):
    'Create context inside directory, returning to original on exit'
    original = os.getcwd()
    os.chdir(directory)
    yield os.getcwd()
    os.chdir(original)

# Close the figure used by plotall() if the variable "output" exists
if "output" in dir() and hasattr(output, 'figure'):
    output.figure.close()
    
output = experiment.List()
rip_size = []
for info in good_mol:
    directory, mol_info, pulls = info[0], info[1:4], info[4]
    construct, conditions, mol = mol_info
    pulls = eval('r_[{}]'.format(pulls))
    
    with ChangeDirectory(path.join(DATA_DIR, directory)):
        try:
            exp = experiment.fromMatch(construct, conditions, mol)
            if len(exp) == 0:
                raise RuntimeError(
                'No pulling files found matching {}\{}\n'.format(
                            directory,'_'.join(mol_info)))
            exp = exp.filter(lambda p: p.info.pull in pulls).has_value(trap_f_atleast=16)
            exp.adjustOffset(to_x=EXTENSION_OFFSET)
            pulls = exp.collapse()
            pulls.fitHandles(x=(HSTART, HEND))
            rip_size.append( pulls.fitRip(RSTART, REND_FORCE)['Lc1'] )
            output += exp
        except Exception as err:
            print err

with ChangeDirectory(ANALYSIS_DIR):
    plt.clf()
    output.plotall(legend=None)
    fname = path.splitext(pulls_filename)[0]
    output.savefig(fname+'_pulls')
    rip_forces = output.findRip(HEND)
    rip_forces = rip_forces[rip_forces[:,1]>MIN_FORCE_CUTOFF]
    rip_forces = rip_forces[rip_forces[:,1]<MAX_FORCE_CUTOFF]

    plt.figure()
    plt.clf()
    plt.hist(rip_forces[:,1], bins=14, range=(MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
    plt.title('Rip force distribution')
    plt.show()
    plt.xlim((MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
    plt.savefig(fname+'_rip_force_histogram')

    plt.figure()
    plt.hist(rip_size)

    print "\nCalculating stats from {} traces".format(len(rip_forces))
    forces = rip_forces[:,1]
    print "Mean rip stats:\n\tforce = {:.2f} ({:.2f})\trip size = {:.2f} ({:.2f})".format(
        mean(forces), std(forces), 
        mean(rip_size), std(rip_size))
    print "Median rip stats:\n\tforce = {:.2f}\tsize = {:.2f}".format(
        median(forces), median(rip_size))