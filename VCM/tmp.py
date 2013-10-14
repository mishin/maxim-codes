# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 14:54:26 2013

@author: Maxim
"""

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
        
if __name__=="__main__":
    test1()