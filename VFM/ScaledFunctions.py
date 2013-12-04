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
        self._histScalingFactor = list()
        self._histX = list()

    def construct_scaling_model(self,x0):
        fH = self.fHigh(x0)
        fL = self.fLow(x0)
        self.x0 = x0
        self.f0 = fH
        scFactor = self._get_scaling_factor(fH,fL)
        _list = zip(self.fHigh._histFpart,self.fLow._histFpart)
        _hist = [self._get_scaling_factor(fh,fl) for fh,fl in _list]
        self._histScalingFactor = np.array(_hist)
        self._histX = self.fHigh._histXpart
        if len(self.fHigh._histFpart)<self._warmup:
            gradSc = self._get_scaling_factor_gradient(x0)
            self.scalingFunc = Taylor1(x0,scFactor,gradSc)
        else:
            self.scalingFunc = RbfMod(self._histX,self._histScalingFactor)

    def get_trust_region_ratio(self,x):
        fSc = self(x)
        fH = self.fHigh(x)
        fSc0 = self.f0
        if fSc==fSc0:
            return float('inf')
        else:
            return (fSc0 - fH) / (fSc0 - fSc)
    
    def initialize_by_points(self,X):
        for x in X[:-1]:
            self.fHigh(x,True)
            self.fLow(x,True)
        self.construct_scaling_model(X[-1])
    
    def get_gradient(self,x):
        fval = self.__call__(x,False)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i] + self._dx
            grad[i] = (self.__call__(X,False)-fval)/self._dx
        return grad


class AdditiveScaling(ScaledFunction):
    def __call__(self,x,save=True):
        return self.fLow(x,save) + self.scalingFunc(x)
    def _get_scaling_factor(self,fH,fL):
        return fH - fL
    def _get_scaling_factor_gradient(self,x):
        fH, gradH = self.fHigh.get_gradient(x,self._dx)
        fL, gradL = self.fLow.get_gradient(x,self._dx)
        return gradH - gradL

class MultiplicativeScaling(ScaledFunction):
    def __call__(self,x,save=True):
        return self.fLow(x,save) * self.scalingFunc(x)
    def _get_scaling_factor(self,fH,fL):
        return fH/fL
    def _get_scaling_factor_gradient(self,x):
        fH, gradH = self.fHigh.get_gradient(x,self._dx)
        fL, gradL = self.fLow.get_gradient(x,self._dx)
        return (gradH*fL - fH*gradL)/(fL*fL)

class HybridScaling:
    def __init__(self,fHigh,fLow,warmUpRun=3,dx=1e-4):
        self._histScalingFactor = None
        self._scalingAdd = AdditiveScaling(fHigh,fLow,warmUpRun,dx)
        self._scalingMult = MultiplicativeScaling(fHigh,fLow,warmUpRun,dx)
        self._w = 0.5

    def __call__(self,x,save=True):
        w = self._w
        return w *self._scalingAdd(x,save) + (1.-w) *self._scalingMult(x,save)

    def construct_scaling_model(self,x0):
        pass
    def get_trust_region_ratio(self):
        pass
    def get_gradient(self):
        pass

class TrustRegionManagement:
    def __init__(self,delta,eta1=0.25,eta2=0.75,eta3=1.25,c1=0.3,c2=2.0):
        self.eta1 = eta1
        self.eta2 = eta2
        self.eta3 = eta3
        self.c1 = c1
        self.c2 = c2
        self.delta0 = delta
    
    def adjust(self,rho,x0,xnew):
        delta0 = self.delta0
        if rho<=self.eta1 or rho>=self.eta3:
            delta = self.c1*delta0
        elif self.eta1<rho<self.eta2:
            delta = delta0
        else:
            delta = delta0 * self._get_gamma(x0,xnew,delta0)
        if rho>0:
            x0new = xnew
        else:
            x0new = x0
        self.delta0 = delta
        return delta, x0new
    
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


def debug1():
    fhigh = forrester
    flow = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5
    fsc = HybridScaling(fhigh,flow,3,1e-6)
    

def run_test1():
    x = np.linspace(0,1,50)
    fhigh = lambda x: forrester(x) + 5
    flow = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5.
    
    fsc = MultiplicativeScaling(fhigh,flow,2,1e-6)
    fsc.construct_scaling_model(0.5)
#    fsc.construct_scaling_model(0.2)
#    fsc.construct_scaling_model(0.8)
#    fsc.construct_scaling_model(0.7)
#    fsc.construct_scaling_model(0.9)
#    fsc.construct_scaling_model(0.5)
    #print fsc.get_trust_region_ratio(0.8)
    print fsc.fHigh._nEval
    #fsc.construct_scaling_model(0.8)
    print fsc.fHigh._nEval

    y1 = np.array([fhigh(_x) for _x in x])
    y2 = np.array([flow(_x) for _x in x])
    y3 = np.array([fsc(_x) for _x in x])
    y4 = np.array([fsc.scalingFunc(_x) for _x in x])
    plt.figure(1)
    plt.hold(True)
    plt.grid(True)
    plt.plot(x,y1,'r-')
    plt.plot(x,y2,'b-')
    plt.plot(x,y3,'g-')
    plt.plot(x,y4,'k-')
    plt.plot(fsc.fHigh._histXpart, fsc.fHigh._histFpart,'ko')
    plt.axis([-0.,1.,-10,20])
    plt.show()

if __name__=="__main__":
    debug1()