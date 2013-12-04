# -*- coding: utf-8 -*-
"""
Created on Wed Dec 04 15:51:15 2013

@author: Maxim
"""

from functionHandling import *

class ScaledFunction:
    """
    fhigh, flow - functions
    fH, fL - function values
    """
    def __init__(self,fHigh,fLow,warmUpRun=3,dx=1e-4):
        self.fHigh = FunctionND(fHigh)
        self.fLow  = FunctionND(fLow)
        self.warmUpRun = warmUpRun
        self._dx = dx
        self._histScaling = list()

    def construct_scaling_model(self,x0):
        fH = self.fHigh(x0,True)
        fL = self.fLow(x0,True)
        self.x0 = x0
        self.fH0 = fH
        scalingFactor = self._get_scaling_factor(fH,fL)
        _list = zip(self.fHigh._histFpart,self.fLow._histFpart)
        _hist = [self._get_scaling_factor(_fH,_fL) for _fH,_fL in _list]
        self._histScaling = np.array([_hist])
        if len(self.fHigh._histFpart)<self.warmUpRun:
            scalingGrad = self._get_scaling_factor_gradient(x0,fH,fL)
            self.scalingFunc = TaylorND(x0,scalingFactor,scalingGrad)
        else:
            self.scalingFunc = RbfMod(self.fHigh._histXpart,self._histScaling)

    def _get_scaling_factor(self,fH,fL):
        if self.scalingType=='add':
            return fH-fL
        elif self.scalingType=='mult':
            return fH/fL
        
    def _get_scaling_factor_gradient(self,x,fH,fL):
        fH, gradH = self.fHigh.get_gradient(x,self._dx,fH)
        fL, gradL = self.fLow.get_gradient(x,self._dx,fL)
        if self.scalingType=='add':
            return gradH - gradL
        elif self.scalingType=='mult':
            return (gradH*fL - gradL*fH) / (fL*fL)

    def get_trust_region_ratio(self,x):
        fSc = self.__call__(x)
        fH = self.fHigh(x,True)
        fSc0 = self.fH0
        if fSc==fSc0:
            return float('inf')
        else:
            return (fSc0 - fH) / (fSc0 - fSc)
    
    def initialize_by_points(self,X):
        for x in X[:-1]:
            self.fHigh(x,True)
            self.fLow(x,True)
        self.construct_scaling_model(X[-1])
    
    def initialize_from_file(self,fileHigh,fileLow=None):
        pass
    
    def get_gradient(self,x):
        fval = self.__call__(x,False)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i] + self._dx
            grad[i] = (self.__call__(X,False)-fval)/self._dx

class AdditiveScaling(ScaledFunction):
    def __call__(self,x,save=True):
        return self.scalingFunc(x) + self.fLow(x,save)

class MultiplicativeScaling(ScaledFunction):
    def __call__(self,x,save=True):
        return self.scalingFunc(x) * self.fLow(x,save)

class HybridScaling:
    def __init__(self):
        pass