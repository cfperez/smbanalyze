# -*- coding: utf-8 -*-
"""
Created on Fri Aug 2, 2013
@author: Christian
"""
import os
import os.path as path
from itertools import groupby
from contextlib import contextmanager
from smbanalyze import experiment, fileIO
from numpy import mean, median, std
import datetime

# DEFAULT: Change the path for your machine!
PROJECT_DIR = r'/Volumes/users2/Force-FRET Project'


@contextmanager
def ChangeDirectory(directory):
    'Create context inside directory, returning to original on exit'
    original = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(original)


def to_date(date_string):
    return datetime.date(*map(int, date_string.split('.')))


class ExperimentBrowser(object):
    DATA_DIR = 'Data'
    ANALYSIS_DIR = 'Analysis'

    def __init__(self, filename=None, exp_type=experiment.Pulling, project_dir=PROJECT_DIR):
        self._exp_type = exp_type
        self.project_dir = project_dir
        self.data_dir = path.join(project_dir,
                                  ExperimentBrowser.DATA_DIR)
        self.analysis_dir = path.join(project_dir,
                                      ExperimentBrowser.ANALYSIS_DIR)
        if filename:
            with open(filename) as fh:
                self._info = [line.strip().split() for line in fh.readlines()]

    def __iter__(self):
        return self._iter()

    def byMolecule(self):
        return self._iter()
        
    def _iter(self, yield_func=None, with_directory=False):
        assert isinstance(with_directory, bool)
        for info in self._info:
            directory, mol_info, pulls = info[0], info[1:4], info[4]
            pulls = eval('r_[{}]'.format(pulls))

            for mol in byMolecule(path.join(self.data_dir, directory),
                                  exp_type=self._exp_type,
                                  fromMatch=mol_info):
                mol = mol.filter(lambda p: p.info.pull in pulls)
                if len(mol) == 0:
                    print 'No pulling files found matching {}\{}\n'.format(
                        directory, '_'.join(mol_info))
                else:
                    if yield_func:
                        yield yield_func(mol, info)
                    else:
                        yield mol

    def with_directory(self):
        return self._iter( yield_func=(lambda m,i: (to_date(i[0]),m)) )

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


def fromTabFile(filename, exp_type=experiment.Pulling, project_dir=PROJECT_DIR):
    return ExperimentBrowser(filename, exp_type=exp_type, project_dir=project_dir)


def groupMolecules(exps):
    group_iter = groupby(exps, lambda p: p.info[:4])
    for info, group in group_iter:
        yield experiment.List(group)

def dir_to_date(date_str):
    try:
        dir_list = map(int, date_str.split('.'))
    except ValueError:
        return None
    dir_list[0] -= 2000
    return tuple(dir_list)

def byMolecule(directory='./', exp_type=experiment.Pulling, fromMatch=('',), **metadata):
    """Return a generator (iterator) which yields experiment.List for all molecules in given directory."""
    assert issubclass(exp_type, experiment.Base)
    assert isinstance(fromMatch, (tuple, list))
    assert isinstance(directory, str)
    with ChangeDirectory(directory):
        flist = filter(exp_type.filenameMatchesType,
                   experiment.filelist(*fromMatch))
    filelist = groupby(flist, lambda f: fileIO.parseFilename(f)[:4])
    for mol, group in filelist:
        with ChangeDirectory(directory):
            dirname = path.basename(directory)
            date_ = dir_to_date(dirname)
            metadata.setdefault('date', date_)
            yield experiment.List(
                map(lambda f: exp_type.fromFile(f, **metadata), group))

def execute(func, filename):
    with ExperimentBrowser(filename) as exp_browser:
        return [func(mol) for mol in exp_browser]
