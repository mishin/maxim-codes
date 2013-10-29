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


def currin(x):
    val1 = 1.0-np.exp(-0.5/x[1])
    val2 = 2300.*x[0]**3 + 1900*x[0]**2 + 2092*x[0]+60
    val3 = 100.*x[0]**3+500*x[0]**2+4*x[0]+20.
    return val1*val2/val3

def currin_low(x):
    x1 = x[0]
    x2 = x[1]
    maxarg = max([0, x2+0.05])
    val1 = currin([x1+1/20, x2+1/20])
    val2 = currin([x1+1/20, maxarg])
    val3 = currin([x1-1/20, x2+1/20])
    val4 = currin([x1-1/20, maxarg])
    return 0.25*(val1+val2+val3+val4)
        
def run_2d_example():
    tol  = 1.0e-4
    lb   = np.array([1.0,1.0])
    ub   = np.array([10.,10.])
    x0   = np.array([5.0,5.0])
    iterMax = 50
    
    xDoe = list()
    xDoe.append([lb[0],lb[1]])
    xDoe.append([ub[0],ub[1]])
    xDoe.append([lb[0],ub[1]])
    xDoe.append([ub[0],lb[1]])
    xDoe.append([(lb[0]+ub[0])/2.,(lb[1]+ub[1])/2.])

    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.3, 2.0)
    err = tol + 1
    niter = 0
    gNotConverged=True
    
    fscaled = ScaledFunction(flow,fhigh,1,'mult')
    gscaled = ScaledFunction(glow,ghigh,1,'mult')
    
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
        rho1,fHiNew = fscaled.get_trust_region_ratio(xnew)
        rho2,gHiNew = gscaled.get_trust_region_ratio(xnew)
        rho = min([rho1,rho2])
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

def currin(x):
    val1 = 1.0-np.exp(-0.5/x[1])
    val2 = 2300.*x[0]**3 + 1900*x[0]**2 + 2092*x[0]+60
    val3 = 100.*x[0]**3+500*x[0]**2+4*x[0]+20.
    return -val1*val2/val3

def currin_low(x):
    x1 = x[0]
    x2 = x[1]
    maxarg = max([0, x2+0.05])
    val1 = currin([x1+1/20, x2+1/20])
    val2 = currin([x1+1/20, maxarg])
    val3 = currin([x1-1/20, x2+1/20])
    val4 = currin([x1-1/20, maxarg])
    return -0.25*(val1+val2+val3+val4)
    
def run_2d_unc():
    tol  = 1.0e-4
    lb   = np.array([0.1,0.1])
    ub   = np.array([1.,1.])
    x0   = np.array([.5,.5])
    iterMax = 50
    
    xDoe = list()
    xDoe.append([lb[0],lb[1]])
    xDoe.append([ub[0],ub[1]])
    xDoe.append([lb[0],ub[1]])
    xDoe.append([ub[0],lb[1]])
    xDoe.append([(lb[0]+ub[0])/2.,(lb[1]+ub[1])/2.])

    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.3, 2.0)
    err = tol + 1
    niter = 0
    
    fscaled = ScaledFunction(currin_low,currin,1,'mult')
    fscaled._initialize_by_doe_points(xDoe)
    
    while err>tol and iterMax>niter:
        fscaled.construct_scaling_model(x0)
        xlb = [max(lb[i],x0[i]-delta) for i in range(len(x0))]
        xub = [min(ub[i],x0[i]+delta) for i in range(len(x0))]
        bnds = [(xlb[0],xub[0]),(xlb[1],xub[1])]
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-6)
        xnew = rslt.x
        fnew = rslt.fun
        rho,fhiNew = fscaled.get_trust_region_ratio(xnew)
        err = np.linalg.norm([x0-xnew])
        delta = trustRegion.adjust(rho,err)
        x0 = xnew
        niter +=1
        print '%.4f\t%.4f\t'%(x0[0],x0[1]),'%.2f\t%.2e\t'%(rho,delta),'%.2e'%err
    fscaled.funcHi.display()
    fscaled.funcLo.display()
    print niter

def exact_solution():
    rslt = minimize(currin,[.5,.5],method='SLSQP',bounds=((0.0,1.),(.0,1.)))
    print rslt

if __name__=="__main__":
    #run_2d_example()
    exact_solution()