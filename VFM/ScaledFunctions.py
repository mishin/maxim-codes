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
        pass

    def construct_scaling_model(self,x0):
        pass

    def _get_scaling_factor(self,fH,fL):
        pass
        
    def _get_scaling_factor_gradient(self,x,fH,fL):
        pass

    def get_trust_region_ratio(self,x):
        pass
    
    def initialize_by_points(self,X):
        pass
    
    def initialize_from_file(self,fileHigh,fileLow=None):
        pass
    
    def get_gradient(self,x):
        pass


class AdditiveScaling(ScaledFunction):
    def __call__(self,x,save=True):
        return self.scalingFunc(x) + self.fLow(x,save)
    def _get_scaling_factor(self,x,fH,fL):
        pass
    def _get_scaling_factor_gradient(self,x,fH,gradH,fL,gradL):
        pass

class MultiplicativeScaling(ScaledFunction):
    def __call__(self,x,save=True):
        return self.scalingFunc(x) * self.fLow(x,save)
    def _get_scaling_factor(self,x,fH,fL):
        pass
    def _get_scaling_factor_gradient(self,x,fH,gradH,fL,gradL):
        pass

class HybridScaling:
    def __init__(self):
        pass
    def __call__(self):
        pass