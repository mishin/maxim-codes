# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 01:20:13 2013

@author: Maxim
"""
import numpy as np
import matplotlib.pyplot as plt

class Function1D:
    def __init__(self,func):
        self.func = func
        self._histX = list()
        self._histF = list()
        self._nEval = 0
        self._nGrad = 0
    
    def __call__(self,x):
        for i,xval in enumerate(self._histX):
            if x==xval:
                return self._histF[i]
        else:
            self._nEval += 1
            f = self.func(x)
            self._histX.append(x)
            self._histF.append(f)
            return f
        
    def get_gradient(self,x,dx=1e-3):
        self._nGrad += 1
        f = self.__call__(x)
        df = self.__call__(x+dx)
        return (df-f)/dx
    
    def get_taylor(self,x,dx=1e-3,fval=None):
        if fval==None:
            fval = self.__call__(x)
        grad = self.get_gradient(x,dx)
        return Taylor1D(x,fval,grad)


class FunctionND:
    def __init__(self,func):
        self.func = func
        self._histX = 0
        self._histF = 0
        self._nEval = 0
        self._nGrad = 0
    
    def __call__(self,x):
        if self._nEval==0:
            self._histX = x
            fval = self.func(x)
            self._histF = fval
        else:
            for i,xval in enumerate(self._histX):
                if all(x==xval):
                    fval = self._histF[i]
                    break
            else:
                self._nEval += 1
                fval = self.func(x)
                self._histX = np.vstack(self._histX,x)
                self._histF = np.hstack([self._histF,fval])
        return fval

    def get_gradient(self,x,dx=1e-3):
        fval = self.__call__(x)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i]+dx
            grad[i] = (self(X,False)-fval)/dx
        return fval, grad
    
    def get_taylor(self,x,dx=1e-3,fval=None):
        if fval==None:
            fval = self.__call__(x)
        grad = self.get_gradient(x,dx)
        return Taylor1D(x,fval,grad)


class Taylor1D:
    def __init__(self,x,f,gradF):
        self.x0 = x
        self.f0 = f
        self.grad = gradF
    def __call__(self,x):
        return self.f0 + (x - self.x0)*self.grad


class TaylorND:
    def __init__(self,x,f,gradF):
        self.x0 = x
        self.f0 = f
        self.grad = gradF
    def __call__(self,x):
        return self.f0 + np.dot((np.array(x)-self.x0),self.grad)


def test1d():
    x = np.linspace(0,1,50)
    func = lambda x: (5.0*x-2.0)**2.0*np.sin(12.*x-4.)
    f = Function1D(func)
    approxF = f.get_taylor(0.9)
    
    plt.figure(1)
    plt.grid(True)
    plt.plot(x,[f(_x) for _x in x],'b-')
    plt.plot(x,[f(_x) for _x in x],'g--')
    plt.hold(True)
    plt.plot(x,[approxF(_x) for _x in x],'r-')
    plt.show()

def test2d():
    func = lambda x: 1

if __name__=="__main__":
    test1d()