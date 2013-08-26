# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:43:38 2013

@author: Christian
"""
import os.path as path
from smbanalyze import experiment, browser
from numpy import r_, mean, median, std, array, ravel
import matplotlib.pyplot as plt


WORKING_DIR = '/Volumes/users2/Force-Fret Project/Analysis/!Rip/'

# curve adjustment parameters
EXTENSION_OFFSET = 1050

# For outliers
MIN_FORCE_CUTOFF = 8.5
MAX_FORCE_CUTOFF = 14

# For FEC fits
HSTART = 900
HEND = 1012
RSTART = 1035
REND_FORCE = 15


def calc_rip_size(pulls):
    p = pulls.collapse()
    p.fitHandles(x=(HSTART, HEND))
    return p.fitRip(RSTART, REND_FORCE)['Lc1']


def calc_rip_sizes_list(pulls):
    pulls.fitHandles(x=(HSTART, HEND))
    pulls.fitRip(RSTART, REND_FORCE)
    return ravel(pulls.get('ripSizes'))


def calc_rip_forces(pulls):
    return pulls.findRip(HEND)[:, 1].tolist()


def process_experiments(pulls):
    'Return rip size and rip forces from list of ExperimentBrowser'
    assert isinstance(pulls, experiment.List)
    print "Calculating {}".format(pulls[0].filename)
    pulls.adjustOffset(to_x=EXTENSION_OFFSET)
    return calc_rip_size(pulls), calc_rip_forces(pulls)


def process(filename, funcs=()):
    """Return loaded experiments and results calculated from funcs"""
    assert isinstance(funcs, tuple)
    all_mol = list(browser.ExperimentBrowser(filename))
    output = [map(func, all_mol) for func in funcs]
    return tuple(all_mol, *output)
    
def plot_distribution(data, title, xlabel):
    plt.figure()
    plt.hist(data, bins=14, normed=True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.show()

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    parser.add_argument("-L", "--load-only", action="store_true")
    args = parser.parse_args()

    # Close the figure used by plotall() if the variable "output" exists
    if "output" in dir() and hasattr(output, 'figure'):
        output.figure.close()

    # output = experiment.List()
    output = []
    if not args.load_only:
        rip_size = []
        rip_forces = []
    filename = path.join(WORKING_DIR, args.file)
    with browser.ExperimentBrowser(filename) as exp_browser:
        for mol in exp_browser:
            mol = mol.has_value(trap_f_atleast=16)
            mol.adjustOffset(to_x=EXTENSION_OFFSET)
            if args.load_only:
                print mol
            else:
                rip_size += [calc_rip_size(mol)]
                # rip_forces += calc_rip_forces(mol)
            output += [mol]

    if not args.load_only:
        rip_size_flat = [arr for arr in flatten(rip_size)]
        with browser.ChangeDirectory(WORKING_DIR):
            plt.clf()
            output.plotall(legend=None)
            figure_name = path.splitext(filename)[0]
            output.savefig(figure_name + '_pulls')

            # plt.figure()
            # plt.clf()
            # plt.hist(rip_forces, bins=14, range=(MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
            # plt.title('Rip force distribution')
            # plt.show()
            # plt.xlim((MIN_FORCE_CUTOFF, MAX_FORCE_CUTOFF))
            # plt.xlabel('Rip force (pN)')
            # plt.ylabel('# of molecules')
            # plt.savefig(figure_name+'_rip_force_histogram')

            plt.figure()
            plt.hist(rip_size)
            plt.title('Rip size distribution')
            plt.xlabel('Rip size (nm)')
            plt.ylabel('# of molecules')
            plt.show()

            print "\nCalculating stats from {} traces".format(len(rip_forces))
            # print "Force = {:.2f} ({:.2f})".format(
                # mean(rip_forces), std(rip_forces))
            print "Rip size = {:.2f} ({:.2f})".format(
                mean(rip_size), std(rip_size))
            # print "Median rip stats:\n\tforce = {:.2f}\tsize = {:.2f}".format(
                # median(rip_forces), median(rip_size))
