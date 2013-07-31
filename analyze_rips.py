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

# For outliers
MIN_FORCE_CUTOFF = 8.5
MAX_FORCE_CUTOFF = 14

# Change the drive letter for your machine!
#DATA_DIR = 'Y:\Force-FRET Project\Data'
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
if "output" in dir():
    output.figure.close()
    
output = experiment.List()
for info in good_mol:
    directory, mol_info, pulls = info[0], info[1:4], info[4]
    construct, conditions, mol = mol_info
    pulls = eval('r_[{}]'.format(pulls))
    
    with ChangeDirectory(path.join(DATA_DIR, directory)):
        try:
            exp = experiment.fromMatch(construct, conditions, mol).is_a(experiment.Pulling)
            if len(exp) == 0:
                raise Exception(
                'No pulling files found matching {}\{}\n'.format(
                            directory,'_'.join(mol_info)))
            exp = exp.filter(lambda p: p.info.pull in pulls)
            exp = exp.filter(lambda p: any(p.trap.f >= 16))
            output += exp
        except Exception as err:
            print err

with ChangeDirectory(ANALYSIS_DIR):
    output.adjustOffset()
    plt.clf()
    output.plotall(legend=None)
    fname = path.splitext(pulls_filename)[0]
    output.savefig(fname+'_pulls')
    rips = output.findRip(980)
    rips = rips[rips[:,1]>MIN_FORCE_CUTOFF]
    rips = rips[rips[:,1]<MAX_FORCE_CUTOFF]

    plt.figure()
    plt.clf()
    plt.hist(rips[:,1], bins=14, range=(MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
    plt.title('Rip force distribution')
    plt.xlim((MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
    plt.savefig(fname+'_rip_force_histogram')

    print "\nCalculating stats from {} rips".format(len(rips))
    print "Mean rip stats: ", mean(rips, axis=0)
    print "Median rip stats: ", median(rips, axis=0)
    print "Std rip stats: ", std(rips, axis=0)
    