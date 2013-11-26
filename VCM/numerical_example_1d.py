# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 21:57:54 2013

@author: Maxim
"""
from scalingFunctions import *

def get_bounds(x,delta,lb,ub):
    bounds = np.zeros([1,2])
    bounds[0,0] = max(lb,x-delta)
    bounds[0,0] = min(ub,x+delta)
    return bounds

def run_example_1d():
    fhigh = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
    flow = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5
    
    tol = 1e-6
    iterMax = 20
    xl = 0.0
    xu = 1.0
    x0 = 0.5
    xdoe = np.array([xl,xu,0.3])
    
    err = tol+1
    nIter = 0
    
    x = np.linspace(-.1,1.1,50)
    yh = np.array(fhigh(_x) for _x in x)
    yl = np.array(flow(_x) for _x in x)
    
    delta = 0.5
    fscaled = ScalingFunction(fhigh,flow,'add',3,True)
    while err>tol and nIter<iterMax:
        fscaled.construct_scaling_model(x0)
        bnds = get_bounds(x0,delta,xl,xu)
        rslt = minimize(fscaled,x0,'SLSQP',bounds=bnds)
        xnew =rslt.x
        fnew = rslt.fun
        rho = fscaled.get_trust_region_ratio(xnew)
        


if __name__=="__main__":
    run_example_1d()