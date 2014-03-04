'''
Progress bar

Created 2/25/14
@author Christian P
'''

import sys
from functools import wraps
from contextlib import contextmanager

DEFAULT_SIZE = 80

def make(status='', size=DEFAULT_SIZE, char='='):
    def progress_bar(p):
        p = min(p,1.0)
        if not progress_bar.done:
            progress_bar.done = (p==1.0)
            prog_amt = '{percent:.0%}' if p<1.0 else 'Done!\n'
            prog_bar = '\r{status} [{bar:{size}s}] ' + prog_amt
            sys.stdout.write(prog_bar.format(
              status=status,bar=char*int(p*size),percent=p,size=size))
            sys.stdout.flush()
    progress_bar.done = False
    return progress_bar

def on_next(N, size=DEFAULT_SIZE, status=''):
    pbar = make(status, size=size)
    N = float(N)
    n = 1
    while True:
        update = (yield pbar(min(n,N)/N))
        n = update*N if update else n+1

def over(iterable, N=None):
    N = N or len(iterable)
    pb = on_next(N) #AutoProgressBar(N)
    for x in iterable:
       next(pb)
       yield x

@contextmanager
def done_on_complete(N, size=DEFAULT_SIZE, status=''):
    pb = on_next(N, size=size, status=status)
    try:
        yield pb
    finally:
        pb.send(True)

def progress_on_call(func, pbar):
    @wraps(func)
    def f_(*args, **kwargs):
        next(pbar)
        return func(*args, **kwargs)
    return f_

#############################
# OOP version of above... 
# Bloated, still useful?
class ProgressBar(object):
    def __init__(self, status='', size=DEFAULT_SIZE, char='='):
        self.status = status
        self.size = size
        self.char = char
        self._done = False
    
    def update(self, p):
        p = min(p,1.0)
        if not self._done:
            self._done = (p==1.0)
            prog_amt = '{percent:.0%}' if p<1.0 else 'Done!\n'
            prog_bar = '\r{status} [{bar:{size}s}] ' + prog_amt
            sys.stdout.write(prog_bar.format(
              status=self.status,bar=self.char*int(p*self.size),percent=p,size=self.size))
            sys.stdout.flush()

    def done(self):
        self.update(1.0)

    def __iter__(self):
        return self
    
    def next(self):
        self.update()

    def __enter__(self):
        return self
    def __exit__(self, type_, value, tb):
        self.done()
        if type_ is not None:
            raise type_(value, tb)

    def __call__(self, iterable):
        size = float(len(iterable))
        for x in enumerate(iterable):
            self.update(x/size)
            yield x

class AutoProgressBar(ProgressBar):
    def __init__(self, maxcalls, status=''):
        super(AutoProgressBar,self).__init__(status=status)
        self._length = maxcalls
        self._count = 0
        
    def update(self, p=None):
        if p is not None:
            self._count = self._length*p
        else:
            if self._count <= self._length:
                self._count += 1
            p = self._count/float(self._length)
        super(AutoProgressBar,self).update(p)

    def __call__(self, iterable):
        for x in iterable:
            self.update()
            yield x
        self.done()