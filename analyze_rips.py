# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:43:38 2013

@author: Christian
"""
import os
from contextlib import contextmanager

DATA_DIR = 'C:\Users\Christian\Documents\FRET\Data'

MOLECULES = '''2013.07.02	SJ2at	s2m3	10pM   1:5,6,7
2013.07.02	SJ2at	s2m4	10pM   1:6,7,8
2013.06.28	SJ2at	s2m1	100pM 1:3
2013.06.14	SJ2at	s5m1	15pM  3:5,6:9'''

good_mol = [m.split() for m in MOLECULES.split('\n')]

os.chdir(DATA_DIR)
print os.getcwd()

@contextmanager
def ChDir(directory):
    os.chdir(directory)
    yield os.getcwd()
    os.chdir('..')
    
output = experiment.List()
for info in good_mol:
    directory, mol_info, pulls = info[0], info[1:4], info[4]
    construct, mol, conditions = mol_info
    pulls = eval('r_[{}]'.format(pulls))
    
    with ChDir(directory):
        try:
            exp = experiment.fromMatch(construct, conditions, mol).is_a(experiment.Pulling)
            exp = exp.filter(lambda p: p.info.pull in pulls)
            exp = exp.filter(lambda p: any(p.trap.f >= 16))
            output += exp
        except Exception as err:
            print err
                
output.adjustOffset()
rips = output.findRip(980)
print "Mean rip stats: ", mean(rips, axis=0)
print "Std rip stats: ", std(rips, axis=0)