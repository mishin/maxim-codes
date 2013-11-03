# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 16:14:35 2013

@author: Maxim
"""

from pyDOE import lhs
from numpy import array, hstack, vstack

def save_samples():
    n = 100
    dim = 2
    path = 'LHC_test_samples_2d.txt'
    xDoe = lhs(dim,n)
    fid = open(path,'wt')
    for xline in xDoe:
        for x in xline:
            fid.write('%.6f\t'%x)
        fid.write('\n')
    fid.close()

def read_samples(path=None):
    if path==None:
        path = 'LHC_test_samples_2d.txt'
    fid = open(path,'rt')
    for i,line in enumerate(fid):
        if line.strip()!='':
            n = len(line)
            xline = array([float(x) for x in line.split()])
            if i==0:
                xDoe=xline
            else:
                xDoe = vstack([xDoe,xline])
    return xDoe

if __name__=="__main__":
    save_samples()
    #print read_samples()