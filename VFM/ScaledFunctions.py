# -*- coding: utf-8 -*-
"""
Created on Wed Dec 04 15:51:15 2013

@author: Maxim
"""

from FunctionHandling import *
import matplotlib.pyplot as plt

class ScaledFunction(object):
    """
    fhigh, flow - functions
    fH, fL - function values
    """
    def __init__(self,fHigh,fLow,warmUpRun=3,dx=1e-4):
        self.fHigh = Function(fHigh)
        self.fLow  = Function(fLow)
        self._warmup = warmUpRun
        self._dx = dx
#        self._histScalingFactor = list()
        self._histX = list()

    def construct_scaling_model(self,x0):
        fH = self.fHigh(x0)
        fL = self.fLow(x0)
        self.x0 = x0
        self.f0 = fH
        
        _list = zip(self.fHigh._histFpart,self.fLow._histFpart)
        _hist = [self._get_scaling_factor(fh,fl) for fh,fl in _list]
        self._histScalingFactor = np.array(_hist)
        self._histX = np.copy(self.fHigh._histXpart)

        if len(self.fHigh._histFpart)<self._warmup:
            scFactor = self._get_scaling_factor(fH,fL)
            gradSc = self._get_scaling_factor_gradient(x0)
            self.scalingFunc = Taylor1(x0,scFactor,gradSc)
        else:
            self.scalingFunc = RbfMod(self._histX,self._histScalingFactor)

    def get_trust_region_ratio(self,x):
        fSc = self(x,True)
        fH = self.fHigh(x,True)
        fSc0 = self.f0
        if fSc==fSc0:
            return float('inf')
        else:
            return (fSc0 - fH) / (fSc0 - fSc)
    
    def initialize_by_points(self,X):
        for x in X:
            self.fHigh(x,True)
            self.fLow(x,True)
    
    def get_gradient(self,x):
        fval = self.__call__(x,False)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i] + self._dx
            grad[i] = (self.__call__(X,False)-fval)/self._dx
        return grad


class AdditiveScaling(ScaledFunction):
    def __call__(self,x,save=False):
        return self.fLow(x,save) + self.scalingFunc(x)

    def _get_scaling_factor(self,fH,fL):
        return fH - fL

    def _get_scaling_factor_gradient(self,x):
        fH, gradH = self.fHigh.get_gradient(x,self._dx)
        fL, gradL = self.fLow.get_gradient(x,self._dx)
        return gradH - gradL

class MultiplicativeScaling(ScaledFunction):
    def __call__(self,x,save=False):
        return self.fLow(x,save) * self.scalingFunc(x)
    def _get_scaling_factor(self,fH,fL):
        return fH/fL
    def _get_scaling_factor_gradient(self,x):
        fH, gradH = self.fHigh.get_gradient(x,self._dx)
        fL, gradL = self.fLow.get_gradient(x,self._dx)
        return (gradH*fL - fH*gradL)/(fL*fL)


class HybridScaling:
    def __init__(self,fHigh,fLow,warmUpRun=3,dx=1e-4,w=0.5):
        self.fHigh = Function(fHigh)
        self.fLow  = Function(fLow)
        self._scalingAdd  = AdditiveScaling(fHigh,fLow,warmUpRun,dx)
        self._scalingMult = MultiplicativeScaling(fHigh,fLow,warmUpRun,dx)
        self._w = w
        self._dx = dx

    def __call__(self,x,save=False):
        w = self._w
        return w *self._scalingAdd(x,save) + (1.0-w) *self._scalingMult(x,save)

    def construct_scaling_model(self,x0):
        self._scalingAdd.construct_scaling_model(x0)
        self._scalingMult.construct_scaling_model(x0)
        self.x0 = self._scalingAdd.x0
        self.f0 = self._scalingAdd.f0

    def get_trust_region_ratio(self,x):
        fSc = self(x,True)
        fH = self._scalingAdd.fHigh(x,True)
        self._scalingMult.fHigh(x,True)
        fSc0 = self.f0
        if fSc==fSc0:
            return float('inf')
        else:
            return (fSc0 - fH) / (fSc0 - fSc)
    
    def initialize_by_points(self,X):
        self._scalingAdd.initialize_by_points(X)
        self._scalingMult.initialize_by_points(X)

    def get_gradient(self):
        fval = self.__call__(x,False)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i] + self._dx
            grad[i] = (self.__call__(X,False)-fval)/self._dx
        return grad

class TrustRegionManagement:
    def __init__(self,delta,eta1=0.25,eta2=0.75,eta3=1.25,c1=0.3,c2=2.0):
        self.eta1 = eta1
        self.eta2 = eta2
        self.eta3 = eta3
        self.c1 = c1
        self.c2 = c2
        self.delta0 = delta
    
    def _adjust_x0(self,rho,x0,xnew):
        if rho>0:
            return xnew
        else:
            return x0
    
    def _adjust_delta(self,rho,x0,xnew):
        delta0 = self.delta0
        if rho<=self.eta1 or rho>=self.eta3:
            return self.c1*delta0
        elif self.eta1<rho<self.eta2:
            return delta0
        else:
            return delta0 * self._get_gamma(x0,xnew,delta0)

    def adjust(self,rho,x0,xnew):
        deltaNew = self._adjust_delta(rho,x0,xnew)
        self.delta0 = deltaNew
        x0new = self._adjust_x0(rho,x0,xnew)
        return deltaNew, x0new
    
    def _get_gamma(self,x0,xnew,delta):
        err = np.linalg.norm(x0-xnew)
        if hasattr(x0,'__iter__'):
            err0 = np.linalg.norm(np.ones(len(x0))*delta)
        else:
            err0 = delta
        if err==err0:
            return self.c2
        else:
            return 1.0


def get_bounds(x,delta,lb,ub):
    if hasattr(x,'__iter__'):
        bnds = np.zeros([len(x),2])
        for i,_x in enumerate(x):
            bnds[i,0] = max([lb[i],_x-delta])
            bnds[i,1] = min([ub[i],_x+delta])
    else:
        bnds = np.zeros([1,2])
        bnds[0,0] = max([lb,x-delta])
        bnds[0,1] = min([ub,x+delta])
    return bnds

def debug1():
    fhigh = forrester
    flow = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5
    fsc = HybridScaling(fhigh,flow,3,1e-6)
    

def run_test1():
    x = np.linspace(0,1,50)
    fhigh = lambda x: forrester(x) + 5
    flow = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5.
    
    fsc = AdditiveScaling(fhigh,flow,0,1e-6)
    fsc.initialize_by_points([0.1,0.6,0.8,0.9])
    fsc.construct_scaling_model(0.45)
    fsc(0.12,False)
    fsc(0.15,False)
    fsc(0.11,False)
    fsc(0.115,False)
    fsc.get_trust_region_ratio(0.11)
    fsc.construct_scaling_model(0.41)

    y1 = np.array([fhigh(_x) for _x in x])
    y2 = np.array([flow(_x) for _x in x])
    y3 = np.array([fsc(_x,False) for _x in x])
#    y4 = np.array([fsc.scalingFunc(_x) for _x in x])
    plt.figure(1)
    plt.hold(True)
    plt.grid(True)
    plt.plot(x,y1,'r-')
    plt.plot(x,y2,'b-')
    plt.plot(x,y3,'g-')
#    plt.plot(x,y4,'k-')
    plt.plot(fsc.fHigh._histXpart, fsc.fHigh._histFpart,'ko')
    plt.axis([-0.,1.,-10,20])
    plt.show()

if __name__=="__main__":
    run_test1()