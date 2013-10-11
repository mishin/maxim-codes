# -*- coding: utf-8 -*-
"""
Created on Tue Oct 08 14:09:50 2013

@author: Maxim
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import Rbf

class TestFunction():
    def __init__(self,func):
        self.func = func
        self.xHist = list()
        self.fHist = list()
        self.nEval = 0
        self.nGrad = 0
    def __call__(self,x):
        self.xHist.append(x)
        f = self.func(x)
        self.fHist.append(f)
        self.nEval += 1
        return f
    def get_gradient(self,x,dx=1e-10):
        self.nGrad += 1
        return 0.5*(self(x+dx) - self(x-dx))/dx
    def display(self,showHist=False):
        print 'Number of function evaluations: %d'%self.nEval
        print 'Number of gradient evaluations: %d'%self.nGrad
        if showHist:
            for xval,fval in zip(self.xHist,self.fHist):
                print xval, '| ',fval


class ScaledFunction():
    def __init__(self,funcLo,funcHi):
        self.funcHi = TestFunction(funcHi)
        self.funcLo = TestFunction(funcLo)
        self.nEval = 0
        self.xPrev = list()
        self.fPrev = list()
        self.betaPrev = list()
        self.taylor = True
    def construct_scaling_model(self,x0):
        self.x0 = x0
        self.xPrev.append(x0)
        if len(self.xPrev)<=3:
            fvalLo, fvalHi, gradLo, gradHi = self._get_data(x0)
            self.beta0 = fvalHi / fvalLo
            self.betaPrev.append(self.beta0)
            self.gradBeta = (gradHi*fvalLo - gradLo*fvalHi) / (fvalLo*fvalLo)
        else:
            self.taylor=False
            fvalLo, fvalHi, gradLo, gradHi = self._get_data(x0, False)
            self.betaPrev.append(fvalHi / fvalLo)
            self.response = Rbf(self.xPrev, self.betaPrev)
        self.fPrev.append(fvalHi)

    def __call__(self,x):
        self.nEval += 1
        if self.taylor:
            beta = self.beta0 + (x-self.x0)*self.gradBeta
            return beta*self.funcLo(x)
        else:
            return self.response(x) * self.funcLo(x)

    def _get_data(self,x,gradient=True):
        fvalLo = self.funcLo(x)
        fvalHi = self.funcHi(x)
        if gradient:
            gradLo = self.funcLo.get_gradient(x)
            gradHi = self.funcHi.get_gradient(x)
        else:
            gradLo = None
            gradHi = None
        return fvalLo, fvalHi, gradLo, gradHi

    def get_thrust_region_ratio(self,x):
        fvalHi = self.funcHi(x)
        fvalScaled = self(x)
        fvalScaled0 = self(self.x0)
        if fvalScaled==fvalScaled0:
            return float('inf')
        else:
            rho = (fvalScaled0 - fvalHi) / (fvalScaled0 - fvalScaled)
            return rho
    
    def display(self):
        print 'Number of function evaluations: %d'%self.nEval

def lofiFunc(x):
    return x+1.0
def hifiFunc(x):
    return x*x - 4.0*x + 2.0

def vcm_test():
    eta1 = 0.25
    eta2 = 0.75
    eta3 = 1.25
    c1 = 0.25
    c2 = 2.0
    tol = 1e-3
    err = tol + 1
    niter = 0
    fscaled = ScaledFunction(lofiFunc,hifiFunc)
    delta = 1.0
    x0 = 0.5
    x = np.linspace(-1,5,50)
    while err>tol:
        fscaled.construct_scaling_model(x0)
        bnds = [(x0-delta, x0+delta)]
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-10)
        xnew = rslt.x
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
        x0 = xnew
        print '%.4f\t%.4f\t%.4f\t%.4f'%(rho, xnew, fnew, delta)
        niter += 1
        plt.figure(1)
        #plt.grid(True)
        plt.hold(True)
        plt.plot(x,fscaled.funcLo(x),'r--')
        plt.plot(x,fscaled.funcHi(x),'b--')
        fsc = fscaled(x)
        plt.plot(x,fsc,'k-')
        plt.plot(fscaled.xPrev,fscaled.fPrev,'ro')
        plt.legend(['Low-fi','Hi-fi','Scaled','Current point','Optimum point'],'lower right')
        plt.axis([-1,5,-20,10])
        plt.show()

    fscaled.funcHi.display()
    fscaled.funcLo.display()

def testFunc1(x):
    return np.sin(x)

def vcm_test2D():
    def _fhigh(x):
        return 4.0*x[0]**2 + x[1]**3 + x[0]*x[1]
    def _ghigh(x):
        return 1.0/x[0] + 1.0/x[1] - 2.0
    def _flow(x):
        return 4.0*(x[0]+0.1)**2 + (x[1]-0.1)**3 + x[0]*x[1]+0.1
    def _glow(x):
        return 1.0/x[0] + 1.0/(x[1]+0.1) - 2.0 - 0.001
    # 0.1 < x[0], x[1] < 10
    fhigh = TestFunction(_fhigh)
    flow = TestFunction(_flow)
    ghigh = TestFunction(_ghigh)
    glow = TestFunction(_glow)


def rbf_test():
    x = np.linspace(-1,15,10)
    y = testFunc1(x)
    rbf = Rbf(x,y)
    xnew = np.linspace(-1,15,100)
    ynew = rbf(xnew)

    plt.figure(1)
    plt.grid(True)
    plt.hold(True)
    plt.plot(x,y,'ro')
    plt.plot(xnew,ynew,'r--')
    plt.plot(xnew,testFunc1(xnew),'k-')
    plt.show()

if __name__=="__main__":
    vcm_test()