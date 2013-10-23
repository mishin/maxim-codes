# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 16:51:31 2013

@author: Maxim
"""

from vcm import *

def vcm_test_1d():
    def _flow(x):
        return x+1.0
    def _fhi(x):
        return x*x - 4.0*x + 2.0
    eta1 = 0.25
    eta2 = 0.75
    eta3 = 1.25
    c1 = 0.3
    c2 = 2.0
    tol = 1.0e-4
    err = tol + 1
    niter = 0
    fscaled = ScaledFunction(_flow,_fhi,3)
    delta = 1.0
    x0 = 1.0
    x = linspace(-1,5,50)
    while err>tol:
        fscaled.construct_scaling_model(x0)
        bnds = [(x0-delta, x0+delta)]
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-10)
        xnew = rslt.x[0]
        fnew = rslt.fun
        rho = fscaled.get_thrust_region_ratio(xnew)
        if rho<=eta1 or rho>=eta3:
            delta *= c1
        elif eta2<rho<eta3:
            if abs(x0-xnew)==delta:
                gamma = c2
            else:
                gamma = 1.0
            delta *= gamma
        err = abs(x0-xnew)
        x0plt = x0
        x0 = xnew
        print 'rho:%.4f\tx:%.4f\tf:%.4f\tdelta:%.4e'%(rho, xnew, fnew, delta)
        niter += 1
    
        plt.figure(1)
        #plt.grid(True)
        plt.hold(True)
        plt.plot(x,fscaled.funcLo(x),'r--')
        plt.plot(x,fscaled.funcHi(x),'b--')
        fsc = [fscaled(xsc) for xsc in x]
        plt.plot(x,fsc,'k-')
        plt.plot(fscaled.xPrev,fscaled.fPrev,'ro')
        plt.plot(array([x0plt,x0plt])+delta,[-50,50],'k--')
        plt.plot(array([x0plt,x0plt])-delta,[-50,50],'k--')
        plt.legend(['Low-fi','Hi-fi','Scaled','Current point','Optimum point'],'lower right')
        plt.axis([-1,5,-20,10])
        plt.show()
        plt.cla()

    fscaled.funcHi.display()
    fscaled.funcLo.display()

if __name__=="__main__":
    vcm_test_1d()