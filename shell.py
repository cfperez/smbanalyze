from matplotlib.cbook import flatten
import matplotlib.pyplot as plt
from numpy import linspace, append, diff
from collections import namedtuple
from smbanalyze import fplot, experiment, curvefit, fec, fcalc
import os

from smbanalyze.experiment import Pulling, on_metadata, group_on
from smbanalyze.fplot import Figure
from smbanalyze.fec import nm_to_nt, Rips

experiment_names = ["Pulling", "experiment", "fplot", "Figure"]
fec_names = ['fec', 'nm_to_nt', 'Rips']
__all__ = experiment_names + fec_names + [ "os", "fcalc",
    "fig", "pretty_rip_sizes", "split_pulls_at_point",
    "pick_pts", "pick_line", "pick_intervals", "Interval",
    "group_on"]

def fig(title_):
    fig_ = fplot.Figure(title_).new()
    fig_.clear()
    plt.title(title_)

Interval = namedtuple('Interval', 'start end')

def pick_pts(num=1):
    return Figure.fromCurrent().pickPoints(num)
    
def pick_intervals(num=1):
    '''
    Return list of Intervals picked from current plot
    '''
    return map(Interval._make, Figure.fromCurrent().pickRegions(num))

def to_date(date_string):
    return datetime.date(*map(int, date_string.split('.')))

def region_to_str(regions):
    return '\n'.join(
        'Region {n}: ({:0.1f}, {:0.1f})'.format(
            n=n, *region) for n,region in enumerate(regions)
        )

def split_pulls_at_point(exps, point):
    '''Return tuple (below,above) distinguished by their relative position above/below "point"''' 
    ext_cutoff, force_cutoff = point
    high, low = experiment.List(), experiment.List()
    for p in exps:
        f_at_ext = p.trap.at(ext=ext_cutoff).f
        if not f_at_ext or f_at_ext > force_cutoff:
            high += [p]
        else:
            low += [p]
    return low, high

def truncate_floats(iterable, places=2):
    fmt_str = '{{:.{}f}}'.format(places)
    return map(lambda s: fmt_str.format(s), iterable)

def pretty_rip_sizes(rip_sizes, helices):
    rips_nm = truncate_floats(rip_sizes, places=1)
    rips_nt = truncate_floats(
        map(nm_to_nt, rip_sizes, helices),
        places=1)

    return ', '.join('{} nt ({} nm)'.format(in_nt, in_nm) 
                    for in_nt,in_nm in zip(rips_nt, rips_nm))
    
def flatten_(arr):
    return list(flatten(arr))

def line_from_two_xy_pts(pts):
    assert len(pts) == 2
    m,b = linear_fit_from_two_xy_pts(pts)
    return lambda x: m*x + b

def linear_fit_from_two_xy_pts(pts):
    assert len(pts) == 2
    pts = sorted(pts)
    x0,y0 = pts[0]
    x1,y1 = pts[1]
    m = (y1-y0)/(x1-x0)
    b = y0 - m*x0
    return m,b
    
def linspace_pts(pts):
    assert len(pts) == 2
    pts = sorted(pts)
    x0,x1 = pts[0][0], pts[1][0]
    return linspace(x0, x1)

def pick_line(plot=True):
    pts = plt.ginput(2)
    ll = line_from_two_xy_pts(pts)
    if plot:
        X = linspace_pts(pts)
        plt.plot(X, ll(X))
    return ll