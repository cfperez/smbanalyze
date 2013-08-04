# -*- coding: utf-8 -*-
"""
Created on Fri Aug 2, 2013
@author: Christian
"""
import os
import os.path as path
import sys
from itertools import groupby
from contextlib import contextmanager
from smbanalyze import experiment, fileIO
from numpy import r_, mean, median, std
import matplotlib.pyplot as plt

# DEFAULT: Change the path for your machine!
PROJECT_DIR = r'/Volumes/users2/Force-FRET Project'

@contextmanager
def ChangeDirectory(directory):
    'Create context inside directory, returning to original on exit'
    original = os.getcwd()
    os.chdir(directory)
    yield os.getcwd()
    os.chdir(original)

class ExperimentBrowser(object):
    DATA_DIR = 'Data'
    ANALYSIS_DIR = 'Analysis'
    def __init__(self, filename, exp_type=experiment.Pulling, project_dir=PROJECT_DIR):
        self._exp_type = exp_type
        self.project_dir = project_dir
        self.data_dir = path.join(project_dir,
                             ExperimentBrowser.DATA_DIR)
        self.analysis_dir = path.join(project_dir,
                             ExperimentBrowser.ANALYSIS_DIR)
        with open(filename) as fh:
            self._info = [line.strip().split() for line in fh.readlines()]

    def __iter__(self):
        for info in self._info:
            directory, mol_info, pulls = info[0], info[1:4], info[4]
            pulls = eval('r_[{}]'.format(pulls))
        
            for mol in byMolecule(path.join(self.data_dir, directory),
                            exp_type=self._exp_type,
                            matching=mol_info):
                mol = mol.filter(lambda p: p.info.pull in pulls)
                if len(mol) == 0:
                    print 'No pulling files found matching {}\{}\n'.format(
                                directory,'_'.join(mol_info))
                else:
                    yield mol

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

def fromTabFile(filename, exp_type=experiment.Pulling, project_dir=PROJECT_DIR):
    return ExperimentBrowser(filename, exp_type=exp_type, project_dir=project_dir)

def byMolecule(directory='.', exp_type=experiment.Pulling, matching=('')):
    with ChangeDirectory(directory):
        flist = filter(exp_type.filenameMatchesType,
            experiment.filelist(*matching))
        flist.sort()
        for mol, group in groupby(flist, lambda f: fileIO.parseFilename(f).mol):
            yield experiment.List(map(exp_type.fromFile, group))


def execute(func, filename):
    with ExperimentBrowser(filename) as exp_browser:
        return [func(mol) for mol in exp_browser]