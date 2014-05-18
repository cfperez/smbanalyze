from smbanalyze import fcalc
from nose.tools import *

def test_fret_counts():
  donor,acceptor = 100,200
  beta,gamma = .13, .9
  for donor,acceptor in zip(range(1,1000,100),range(1,1000,100)):
    yield check_fret_counts, donor, acceptor, beta, gamma

def check_fret_counts(donor,acceptor,beta,gamma):
  a_beta = acceptor-donor*beta
  d_,a_,fret=fcalc.fret_counts(donor, acceptor, beta, gamma)
  eq_(donor, d_)
  eq_(a_beta, a_)
  eq_(a_beta/(a_beta+donor*gamma), fret)