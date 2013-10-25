# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 23:08:07 2013

@author: Maxim
"""

from scipy.optimize import minimize
import numpy as np

def test1():
    def obj(x):
        return x[0]+x[1]
    
    def cnstr1(x):
        return x[0]-1.0
    
    def cnstr2(x):
        return x[1]-2.0
    
    bnds = ((-20,20),(-20,20))
    cnstr = ({'type':'ineq','fun':cnstr1},{'type':'ineq','fun':cnstr2})
    print minimize(obj,[5.0,5.0],method='SLSQP',constraints=cnstr,bounds=bnds)
    
if __name__=="__main__":
    test1()