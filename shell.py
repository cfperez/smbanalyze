from matplotlib.cbook import flatten
import matplotlib.pyplot as plt
from numpy import linspace, append, diff
from collections import namedtuple
from smbanalyze import fplot, experiment, curvefit, fec, fcalc, db
import os
import datetime

from smbanalyze.experiment import Pulling, on_metadata, group_dict
from smbanalyze.fplot import Figure
from smbanalyze.fec import nm_to_nt, Rips

_modules_ = ["fplot", "fileIO", "experiment", "curvefit", "fec", "fcalc", "db"]
_experiment_names_ = ["Pulling", "Figure"]
fec_names = ['fec', 'nm_to_nt', 'Rips']
__all__ = _modules_ + _experiment_names_ + fec_names \
    + [ "os", "fcalc", "fig", "pretty_rip_sizes", 
    "split_pulls_at_point", "pick_pts", "pick_line",
    "pick_intervals", "Interval", "group_dict", "to_date",
    "savefig", "plot_segmented", "reload_all", "today"]

Interval = namedtuple('Interval', 'start end')

def today():
    return datetime.datetime.today()

def reload_all():
    reload(fplot)
    reload(fileIO)
    reload(experiment)
    reload(fec)
    reload(db)
    reload(shell)

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
    ratio = p.metadata['sampling_ratio']
    exp_time = (p.metadata['fret.exposurems']/1000.)
    hold(True)
    subplot(211)
    plt.title(title)
    plot(fret.time, fret.fret, 'k:')
    subplot(212)
    for _ in exp:
        plot(_.trap.ext, _.trap.f, 'k:', linewidth=1)
    for start_n in range(len(fret)):
        subplot(211)
        time,donor,acc,E = fret[start_n:start_n+2]
        plot(time,E,'o',markersize=8)
        #plt.autoscale(tight=True)
        subplot(212)
        x,f,sep = trap[start_n*ratio:(start_n+1)*ratio+1].T
        plot(x,f,'-', linewidth=3)
    subplot(211)
    xlim(0,time+exp_time)
    ylim(-0.05,1.1)
    xlabel('Time (s)')
    ylabel('FRET')
    subplot(212)
    xlabel('Extension (nm)')
    ylabel('Force (pN)')
    plt.autoscale(tight=True)

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

def transposed(iterable):
    return [[y[n] for y in iterable] for n in range(len(iterable[0]))]

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