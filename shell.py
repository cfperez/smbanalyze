from matplotlib.cbook import flatten
import matplotlib.pyplot as plt
from numpy import linspace, append, diff
from collections import namedtuple
from smbanalyze import fplot, experiment, curvefit


from smbanalyze.experiment import Pulling
from smbanalyze.fplot import Figure
from smbanalyze.fec import nm_to_nt, Rips

experiment_names = ["Pulling"]
fec_names = ['nm_to_nt', 'Rips']
__all__ = experiment_names + [
    "fig", "nm_to_nt", "pretty_rip_sizes", "split_pulls_at_point",
    "pick_line", "pick_interval", "Interval"]

def fig(title_):
    fig_ = fplot.Figure(title_).new()
    fig_.clear()
    plt.title(title_)

Interval = namedtuple('Interval', 'start end')

def pick_intervals(num=1):
    '''
    Return list of Intervals picked from current plot
    '''
    return map(Interval._make, Figure.fromCurrent().pickRegions(num))

def lc_from_region_fit(fit):
    return [fit[param] for param in fit if param.startswith('Lc')]

def rip_sizes_from_region_fit(fit):
    rips = lc_from_region_fit(fit)[1:]
    return append(rips[0], diff(rips))
    
def region_to_str(regions):
    return '\n'.join(
        'Region {n}: ({:0.1f}, {:0.1f})'.format(
            n=n, *region) for n,region in enumerate(regions)
        )

def split_pulls_at_point(exps, point):
    '''Return tuple of (low,high) of experiments "exps" split at point (ext,force)''' 
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
    
    return "In nm: " + ', '.join(truncate_floats(rip_sizes)) + "\n" + \
         "In nt: " + ', '.join(truncate_floats(map(nm_to_nt, rip_sizes, helices)))

def plot_fit_attr_all(all_mol, attr, bins=10):
    attr_out = [get_fit_attr_from_mol(mol, attr) for mol in all_mol]
    # attr_all = list(flatten(attr_out))
    plt.figure(attr)
    # clf()
    # hist(attr_all, bins=bins, range=(15,50))
    return attr_out

def get_fit_attr_from_mol(mol, attr, bins=10):
    return [p.rips[0][attr] for p in mol]

def flatten_(arr):
    return list(flatten(arr))

def lower_limit_tester(handle_limits, upper_limit):
    def tester(pull, **fitOptions):
        output = []
        for handle in handle_limits:
            pull.fitHandles(handle)
            output.append(pull.fitRip(*upper_limit, **fitOptions))
        return output
    return tester

def region_tester(handle_limits, upper_limit, **fitOptions):
    def tester(pull, **fitOptions):
        output = []
        for handle in handle_limits:
            trap = pull.trap
            handle_mask = trap.maskFromLimits(handle, None)
            upper_mask = trap.maskFromLimits(*upper_limit)
            output.append(
                curvefit.fitWLC_masks(trap.ext, trap.f, [handle_mask, upper_mask], **fitOptions)
            )
        return output
    return tester

def handle_limits(lower_range, upper):
    return [ (lower, upper) for lower in range(*lower_range)]

def fitBoth(pull, hlimit, ulimit):
    pull.fitHandles(hlimit)
    old_fit = pull.fitRip(*ulimit)
    trap = pull.trap
    handles = trap.select(hlimit)
    upper = trap.select(*ulimit)
    new_fit = curvefit.fitWLCrip(handles.ext, handles.f, upper.ext, upper.f, Lp1=1.2)
    return old_fit, new_fit


def fitMultipleLimits(pull, handle_limits, upper_limits):
    for handle, upper in zip(handle_limits, upper_limits):
        old, new = fitBoth(pull, handle, upper)
        yield old, new

def fitCheck(pull, lower_handle_range, upper_handle):
    hlimits = [(lower, upper_handle) for lower in lower_handle_range]
    ulimits = [(1000,15)]*len(hlimits)
    return [[old,new] for old, new in fitMultipleLimits(pull, hlimits, ulimits)]

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
    
def linspace_from_two_xy_pts(pts):
    assert len(pts) == 2
    pts = sorted(pts)
    x0,x1 = pts[0][0], pts[1][0]
    return linspace(x0, x1)

def pick_line():
    pts = plt.ginput(2)
    ll = line_from_two_xy_pts(pts)
    return ll