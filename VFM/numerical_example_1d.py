# -*- coding: utf-8 -*-
"""
Created on Thu Dec 05 20:49:28 2013

@author: Maxim
"""

from ScaledFunctions import *
from scipy.optimize import minimize

fhigh = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
flow  = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5.

def run_test_additive():
    xL, x0, xU = 0.0, 0.5, 1.0
    delta = 0.5
    xDOE = [0.0,1.0]
    fsc = AdditiveScaling(fhigh,flow,0,1e-6)
    fsc.initialize_by_points(xDOE)
    trReg = TrustRegionManagement(delta)
    tol = 1e-3
    err = tol+1
    iterMax = 50
    nIter = 0
    x = np.linspace(0,1,100)
    yh = np.array([fhigh(_x) for _x in x])
    yl = np.array([flow(_x) for _x in x])
    while err>tol and nIter<iterMax:
        fsc.construct_scaling_model(x0)
        

        
        bnds = get_bounds(x0,delta,xL,xU)
        rslt = minimize(lambda _x: fsc(_x,False),x0,method='SLSQP',bounds=bnds)
        xnew = rslt.x[0]
        fnew = rslt.fun
        rho = fsc.get_trust_region_ratio(xnew)
        err = abs(x0-xnew)
        delta,x0 = trReg.adjust(rho,x0,xnew)
        nIter += 1
        print '%.4f\t%.4f\t%.2e\t%.4f\t%.4f'%(xnew,x0, err, fnew, rho)
    plt.figure(1)
    plt.hold(True)
    plt.plot(x,yh,'b')
    plt.plot(x,yl,'r')
    plt.plot(x, [fsc(_x) for _x in x],'k')
    plt.plot(x0,fsc.fHigh(x0,False),'ro')
    plt.plot(fsc._histX,fsc.fHigh._histFpart,'go')
    plt.axis([0,1,-10,20])
    plt.show()
    plt.cla()

if __name__=="__main__":
    run_test_additive()