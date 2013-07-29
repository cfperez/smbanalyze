# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:43:38 2013

@author: Christian
"""
import os
from contextlib import contextmanager

# For outliers
MIN_FORCE_CUTOFF = 8.5
MAX_FORCE_CUTOFF = 15

# Change the drive letter for your machine!
DATA_DIR = 'Y:\Force-FRET Project\Data'
os.chdir(DATA_DIR)
print os.getcwd()

MOLECULES = '''2013.07.02	SJ2at 10pM    s2m3   1:5,6,7
2013.07.02	SJ2at    10pM    s2m4    1:6,7,8
2013.06.28	SJ2at    100pM  s2m1    1:3
2013.06.14	SJ2at    15pM    s5m1    3:5,6:9
2013.04.19  SJF2at  200pM   s1m2    1:5
2013.04.19  SJF2at  200pM   s1m6    1:7
2013.04.19  SJF2at  100pM   s1m1    1:4
2013.04.19  SJF2at  100pM   s1m2    1:4
2013.04.18  SJ6v2at 200pM   s1m2    1:6
2013.04.18  SJ6v2at 200pM   s1m3    1:5
2013.04.18  SJ6v2at 200pM   s1m4    1:6
2013.05.24  SJ2at   20pM    s1m1    1:8
2013.05.24  SJ2at   20pM    s2m1    1:5
2013.05.31  SJ2at  20pM    s4m1    1:13
2013.05.31  SJ2at  100pM   s1m2    1:6
2013.05.31  SJ2at  100pM   s5m1    1:7
2013.05.31  SJ2at  100pM   s5m3    1:5
2013.06.05  SJ2at  100pM   s1m2    1:8
2013.06.05  SJ2at  100pM   s5m1    1:6
2013.06.06  SJ2at  100pM   s2m4    1:10
2013.06.06  SJ2at  50pM   s3m3    1:10,11
2013.06.06  SJ2at  100pM   s5m2    1:6
2013.06.06  SJ2at  50pM   s3m1    1:5,6
2013.06.06  SJ2at  50pM   s3m2    1:11'''

good_mol = [m.split() for m in MOLECULES.split('\n')]


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
    
    with ChangeDirectory(directory):
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

output.adjustOffset()
output.plotall(legend=None)
rips = output.findRip(980)
rips = rips[rips[:,1]>MIN_FORCE_CUTOFF]
rips = rips[rips[:,1]<MAX_FORCE_CUTOFF]

figure()
pyplot.hist(rips[:,1])
title('Rip force distribution')

print "\nCalculating stats from {} rips".format(len(rips))
print "Mean rip stats: ", mean(rips, axis=0)
print "Median rip stats: ", median(rips, axis=0)
print "Std rip stats: ", std(rips, axis=0)
