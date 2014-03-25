import experiment
from shell import split_pulls_at_point, transposed
from scipy.stats import norm, beta
from numpy import sqrt, mean, where
from matplotlib.pyplot import figure,errorbar,ylim,xlim

FILENAME_TOKEN = 'refold'
REFOLD_FILENAME_INFO = ('refold_time', 'series')

def by_time(exps):
    return experiment.group_by(exps, 
        experiment.on_metadata('trap.refolding_time'))

def pretty_str(time_dict):
    outstr = ''
    for key,val in time_dict.iteritems():
        outstr += '\n'.join(
         ("{}s refolding | N = {}".format(key, len(val)),
        '-'*25,
        str(val),
        '\n'))
    return outstr

def find(condition):
    return where(condition)[0]

def avg_time_until_f(exps, force):
    return mean([find(p.trap.f <8)[-1] * p.metadata['trap.sampling_time'] for p in exps])

def adjust_refold_time(by_times, until_force):
    return {t+avg_time_until_f(exps, until_force):exps for t,exps in by_times.items()}

def count_bound(by_times, split_pt):
    return {t:map(len, split_pulls_at_point(exps,split_pt)) for t,exps in by_times.items()}

def wilson_score_z(z=None, confidence=0.683):
    z = z or norm.ppf(0.5+confidence/2)
    def wilson_score(k, n):
        p = float(k)/n
        '''Returns binomial error estimate (lower, upper) given fraction p and samples n'''
        score = lambda sigma: 1/(1 + z**2/n) * (p + z**2/(2*n) + sigma)
        sigma = z*sqrt( p*(1-p)/n + z**2/4/n**2)
        return score(-sigma), score(sigma)
    return wilson_score

def binomial_error(k, n, confidence=0.683):
    alpha = 1-confidence
    return 1-beta.isf(alpha/2.,n-k,k+1), 1-beta.isf(1-alpha/2.,n-k+1,k)

# Using default of 1 standard deviation
wilson_score = wilson_score_z()

def calc_bound_prob(num_bound):
    rtimes,prob,errors = [],[],[]
    for rtime,(unbound,bound) in num_bound.items():
        total = float(unbound+bound)
        frac = bound/float(total)
        l,u = wilson_score(bound,total)
        rtimes.append(rtime)
        prob.append(frac)
        errors.append(map(abs, [frac-l, u-frac]))
    return rtimes, prob, errors

def plot(rtimes, prob, errors):
    errorbar(rtimes, prob, ecolor='b', fmt='o', yerr=transposed(errors))
    ylim(-.05,1.05)
    xlim(0.1,40)
