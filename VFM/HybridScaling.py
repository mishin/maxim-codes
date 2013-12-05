# -*- coding: utf-8 -*-
"""
Created on Thu Dec 05 22:45:25 2013

@author: Maxim
"""
from ScaledFunctions import *


class HybridScaling:
    pass


class _HybridScaling:
    def __init__(self,fHigh,fLow,warmUpRun=3,dx=1e-4,w=0.5):
        self._histScalingFactor = None
        self._scalingAdd  = AdditiveScaling(None,None,warmUpRun,dx)
        self._scalingMult = MultiplicativeScaling(None,None,warmUpRun,dx)
        self.fHigh = Function(fHigh)
        self.fLow  = Function(fLow)
        self._scalingAdd.fHigh = self.fHigh
        self._scalingAdd.fLow  = self.fLow
        self._scalingMult.fHigh = self.fHigh
        self._scalingMult.fLow = self.fLow
        self._w = w
        self._dx = dx

    def __call__(self,x,save=True):
        w = self._w
        return w *self._scalingAdd(x,save) + (1.-w) *self._scalingMult(x,save)


    def construct_scaling_model(self,x0):
        self._scalingAdd.construct_scaling_model(x0)
        self._scalingMult.construct_scaling_model(x0)
        self.f0 = self._scalingAdd.f0

    def get_trust_region_ratio(self,x):
        fSc = self(x,True)
        fH = self._call_high(x)
        fSc0 = self.f0
        if fSc==fSc0:
            return float('inf')
        else:
            return (fSc0 - fH) / (fSc0 - fSc)
    
    def initialize_by_points(self,X):
        for x in X:
            self._scalingAdd.fHigh(x,True)
            self._scalingAdd.fLow(x,True)
        self._copy_history()

    def get_gradient(self):
        fval = self.__call__(x,False)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i] + self._dx
            grad[i] = (self.__call__(X,False)-fval)/self._dx
        return grad