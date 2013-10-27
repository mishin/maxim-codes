# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 18:03:14 2013

@author: Maxim
"""
import numpy as np
from vcm import *

"""
exact solution: x1=0.8842, x2=1.1507, f=5.6683
"""

def flow(x):
    return 4.*(x[0]+0.1)**2.+(x[1]-0.1)**3.+x[0]*x[1]+0.1
def fhigh(x):
    return 4.*x[0]**2.+x[1]**3.+x[0]*x[1]
def glow(x):
    return -(1./x[0]+1./(x[1]+0.1)-2.001)
def ghigh(x):
    return -(1./x[0]+1./x[1]-2.)

def run_2d_example():
    tol  = 1.0e-4
    lb   = np.array([0.1,0.1])
    ub   = np.array([10.,10.])
    x0   = np.array([1.5,1.5])
    iterMax = 50
    
    xDoe = list()
    xDoe.append([lb[0],lb[1]])
    xDoe.append([ub[0],ub[1]])
    xDoe.append([lb[0],ub[1]])
#    xDoe.append([ub[0],lb[1]])
#    xDoe.append([(lb[0]+ub[0])/2.,(lb[1]+ub[1])/2.])

    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.3, 2.0)
    err = tol + 1
    niter = 0
    gNotConverged=True
    
    fscaled = ScaledFunction(flow,fhigh,0,'add')
    gscaled = ScaledFunction(glow,ghigh,0,'add')
    
    fscaled._initialize_by_doe_points(xDoe)
    gscaled._initialize_by_doe_points(xDoe)
    
    while err>tol and iterMax>niter and gNotConverged:
        fscaled.construct_scaling_model(x0)
        gscaled.construct_scaling_model(x0)
        bnds = [(x0[0]-delta,x0[0]+delta),(x0[1]-delta,x0[1]+delta)]
        cnstr = {'type':'ineq','fun':gscaled}
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-6,
                        constraints=cnstr)
        xnew = rslt.x
        fnew = rslt.fun
        rho = min(fscaled.get_thrust_region_ratio(xnew),
                  gscaled.get_thrust_region_ratio(xnew))
        err = np.linalg.norm([x0-xnew])
        delta = trustRegion.adjust(rho,err)
        x0 = xnew
        niter +=1
        print '%.4f\t%.4f\t'%(x0[0],x0[1]),'%.2f\t%.2e\t'%(rho,delta),'%.2e'%err,'%.4f\t%.4f\t'%(glow(xnew),ghigh(xnew))
    gscaled.funcHi.display()
    fscaled.funcHi.display()
    gscaled.funcLo.display()
    fscaled.funcLo.display()
    print niter
if __name__=="__main__":
    run_2d_example()