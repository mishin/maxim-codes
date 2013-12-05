# -*- coding: utf-8 -*-
"""
Created on Fri Dec 06 00:14:47 2013

@author: Maxim
"""
import numpy as np
from scipy.interpolate import Rbf

class Function(object):
    def __init__(self,func,dx=1e-6):
        self.func = func
        self._dx = dx
        self._nEval = 0
        self._nGrad = 0
    
    def __call__(self,x):
        self._nEval += 1
        return self.func(x)
    
    def get_gradient(self,x,fval=None,dx=None):
        self._nGrad += 1
        if dx==None:
            dx = self._dx
        if fval==None:
            fval = self(x)
        if not hasattr(x,'__iter__'):
            x = np.array([x])
            return self._get_gradient(x,fval,dx)[0]
        else:
            return self._get_gradient(x,fval,dx)
    
    def _get_gradient(self,x,fval,dx):
            grad = np.zeros(len(x))
            for i in range(len(x)):
                X = np.copy(x)
                X[i] = X[i]+dx
                grad[i] = (self(X)-fval)/dx
            return fval, grad


class Taylor1:
    def __init__(self,x,f,gradF):
        self.x0 = x
        self.f0 = f
        self.grad = gradF
    def __call__(self,x):
        return self.f0 + np.dot((np.array(x)-self.x0),self.grad)


class RbfMod():
    def __init__(self,x,y):
        if x.ndim>1:
            x = np.transpose(x)
            x = self._get_tuple(x)
            args = x + (y,)
        else:
            args = (x,y)
        self.rbf = Rbf(*args)

    def __call__(self,x):
        if hasattr(x,'__iter__'):
            x = self._get_tuple(x)
            return self.rbf(*x)
        else:
            return self.rbf(x)

    def _get_tuple(self,xArray):
        xTuple = tuple()
        for x in xArray:
            xTuple += (x,)
        return xTuple