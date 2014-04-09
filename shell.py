from matplotlib.cbook import flatten
import matplotlib.pyplot as plt
from numpy import linspace, append, diff
from collections import namedtuple
import fplot, experiment, curvefit, fec, fcalc, db, fileIO
import os
import datetime

_modules_ = ["fplot", "fileIO", "experiment", "curvefit",
         "fec", "fcalc", "db", "os"]
# from importlib import import_module
# _imported_ = map(import_module, _modules_)

# for modname in _modules_:
    # mod = import_module(modname)
    # _all_ = mod.__all__
from fplot import Figure
from experiment import *
from fec import nm_to_nt, Rips
from date import today, date, to_date

_names_ = experiment.__all__

_date_ = ['today', 'date', 'to_date']

fec_names = ['fec', 'nm_to_nt', 'Rips']

__all__ = _modules_ + _names_ + fec_names + _date_ \
    + [ "fig", "pretty_rip_sizes", "pick_pts", "pick_pt", "pick_line",
    "pick_intervals", "Interval", "group_by", "to_date", "transposed",
    "savefig", "plot_segmented", "reload_all"]

Interval = namedtuple('Interval', 'start end')

def reload_all():
    reload(fplot)
    reload(fileIO)
    reload(experiment)
    reload(fec)
    reload(db)

def group_by(iterable, keyfunc):
    return {key: List(p) for key,p in groupby(iterable, keyfunc)}

def savefig(fname=None, **kwargs):
    fname = fname or plt.gca().get_title()
    plt.savefig(fname, transparent=True, **kwargs)

def fig(title_):
    fig_ = fplot.Figure(title_).new()
    fig_.clear()
    plt.title(title_)
    return fig_

def hspan(start, stop, color=None, alpha=.2):
    color = color or next(fplot.COLOR)
    plt.axhspan(start, stop, facecolor=color, alpha=alpha)

def plot_segmented(p, title='', exp=[]):
    trap,fret = p.trap,p.fret
    ratio = p['sampling_ratio']
    exp_time = p['fret.exposurems']/1000.
    plt.hold(True)
    plt.subplot(211)
    plt.title(title)
    plt.plot(fret.time, fret.fret, 'k:')
    plt.subplot(212)
    for _ in exp:
        plt.plot(_.trap.ext, _.trap.f, 'k:', linewidth=1)
    for start_n in range(len(fret)):
        plt.subplot(211)
        time,donor,acc,E = fret[start_n:start_n+2]
        plt.plot(time, E, 'o', markersize=8)
        plt.subplot(212)
        x,f,sep = trap[start_n*ratio:(start_n+1)*ratio+1].T
        plt.plot(x, f, '-', linewidth=3)
    plt.subplot(211)
    plt.xlim(0,time+exp_time)
    plt.ylim(-0.05,1.1)
    plt.xlabel('Time (s)')
    plt.ylabel('FRET')
    plt.subplot(212)
    plt.xlabel('Extension (nm)')
    plt.ylabel('Force (pN)')
    plt.autoscale(tight=True)

def pick_pt():
    return pick_pts(1)[0]

def pick_pts(num=1):
    return Figure.fromCurrent().pickPoints(num)
    
def pick_intervals(num=1):
    '''
    Return list of Intervals picked from current plot
    '''
    return map(Interval._make, Figure.fromCurrent().pickRegions(num))

def region_to_str(regions):
    return '\n'.join(
        'Region {n}: ({:0.1f}, {:0.1f})'.format(
            n=n, *region) for n,region in enumerate(regions)
        )

def transposed(iterable):
    return [[y[n] for y in iterable] for n in range(len(iterable[0]))]

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