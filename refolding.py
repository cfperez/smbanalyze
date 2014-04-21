import experiment
#from explist import List
from scipy.stats import norm, beta
from numpy import sqrt, mean, where, isnan
from matplotlib.pyplot import figure,errorbar,ylim,xlim,gca

FILENAME_TOKEN = 'refold'
REFOLD_FILENAME_INFO = ('refold_time', 'series')

def by_time(exps, key='trap.refolding_time'):
    keyfunc = experiment.on_metadata(key)
    return experiment.group_by(sorted(exps, key=keyfunc), keyfunc)

def prettyprint(time_dict):
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
    return mean([find(p.trap.f<force)[-1] * p.metadata['trap.sampling_time'] for p in exps])

def adjust_refold_time(by_times, until_force):
    return {t+avg_time_until_f(exps, until_force):exps for t,exps in by_times.items()}

def count_bind_tpp(exps):
    return sum(p['bind_tpp']>0 for p in exps), len(exps)

def count_bind_strong(exps):
    return sum(p['bind_tpp']==2 for p in exps), len(exps)

def count(by_times, using=count_bind_tpp):
    return {t:using(exps) for t,exps in by_times.items()}

def count_bound(by_times, split_pt):
    return {t:map(len, experiment.split_pulls_at_point(exps,split_pt)) for t,exps in by_times.items()}

def wilson_score_z(z=0, confidence=0.683):
    z = float(z) or norm.ppf(0.5+confidence/2)
    def wilson_score(k, n):
        p = float(k)/n
        '''Returns binomial error estimate (lower, upper) given fraction p and samples n'''
        score = lambda sigma: 1/(1 + z**2/n) * (p + z**2/(2*n) + sigma)
        sigma = z*sqrt( p*(1-p)/n + z**2/4/n**2)
        return score(-sigma), score(sigma)
    return wilson_score

def binomial_error(k, n, confidence=0.683):
    '''Return (lower,upper) confidence interval of ratio k/n using
    Clopper-Pearson interval (Beta distribution)
    '''
    alpha = 1-confidence
    return beta.ppf(alpha/2.,k,n-k+1), beta.ppf(1-alpha/2.,k+1,n-k)

# Using default of 1 standard deviation
wilson_score = wilson_score_z(1.0)

def calc_bound_prob(num_bound, confidence=.683, error_func=binomial_error):
    rtimes,prob,errors = [],[],[]
    for rtime,(bound,total) in num_bound.items():
        frac = bound/float(total)
        l,u = error_func(bound, total, confidence=confidence)
        if isnan(l):
            l = 0.
        if isnan(u):
            u = 1.
        rtimes.append(rtime)
        prob.append(frac)
        errors.append(map(abs, [frac-l, u-frac]))
    return rtimes, prob, errors

def transposed(iterable):
    return [[y[n] for y in iterable] for n in range(len(iterable[0]))]

def plot(rtimes, prob, errors, color='b'):
    errorbar(rtimes, prob, ecolor=color, fmt='o', yerr=transposed(errors))
    gca().set_xscale('log')
    ylim(-.05,1.05)
    xlim(0.1,40)
