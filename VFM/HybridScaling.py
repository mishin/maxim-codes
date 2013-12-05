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
        self._scalingAdd = AdditiveScaling(fHigh,fLow,warmUpRun,dx)
        self._scalingMult = MultiplicativeScaling(fHigh,fLow,warmUpRun,dx)
        self._w = w

    def __call__(self,x,save=True):
        w = self._w
        return w *self._scalingAdd(x,save) + (1.-w) *self._scalingMult(x,save)
    
    def _call_high(self,x,save=True):
        fH = self._scalingAdd.fHigh(x,save)
        self._copy_history()
        return fH

    def construct_scaling_model(self,x0):
        self._scalingAdd.construct_scaling_model(x0)
        self._copy_history()
        self._scalingMult.construct_scaling_model(x0)
        self.f0 = self._scalingAdd.f0
    
    def _copy_history(self):
        self._scalingMult.fHigh._histFfull = self._scalingAdd.fHigh._histFfull
        self._scalingMult.fHigh._histFpart = self._scalingAdd.fHigh._histFpart
        self._scalingMult.fHigh._histXfull = self._scalingAdd.fHigh._histXfull
        self._scalingMult.fHigh._histXpart = self._scalingAdd.fHigh._histXpart
        self._scalingMult.fLow._histFfull = self._scalingAdd.fLow._histFfull
        self._scalingMult.fLow._histFpart = self._scalingAdd.fLow._histFpart
        self._scalingMult.fLow._histXfull = self._scalingAdd.fLow._histXfull
        self._scalingMult.fLow._histXpart = self._scalingAdd.fLow._histXpart
        self._scalingMult.fHigh._nEval = self._scalingAdd.fHigh._nEval
        self._scalingMult.fLow._nEval = self._scalingAdd.fLow._nEval

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