# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 14:54:26 2013

@author: Maxim
"""
import numpy as np

def test1():
    import numpy as np
    from scipy.interpolate import Rbf
    
    a = np.arange(0,16)
    a = a.reshape(4,4)
    x = tuple()
    for line in a:
        x = x + (line,)
    x = x + (np.array([3,4,5,1.]),)
    c,d,e,f,g = x
    r = Rbf(*x)

def test2():
    x = np.array([1,2,3,4,5,6,7])
    n = len(x)
    if n%2==0:
        Au = x[:n/2]
        Al = x[n/2:]
    else:
        Au = x[:int(n/2)+1]
        Al = np.hstack([x[0],x[int(n/2)+1:]])
    print Au
    print Al

if __name__=="__main__":
    test2()