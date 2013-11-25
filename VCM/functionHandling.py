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
        self._histXfull = list()
        self._histFfull = list()
        self._histXpart = list()
        self._histFpart = list()
        self._nEval = 0
        self._nGrad = 0
    
    def __call__(self,x,save=True):
        for i,xval in enumerate(self._histXfull):
            if x==xval:
                return self._histFfull[i]
        else:
            self._nEval += 1
            f = self.func(x)
            self._histXfull.append(x)
            self._histFfull.append(f)
            if save:
                self._histXpart.append(x)
                self._histFpart.append(x)
            return f
        
    def get_gradient(self,x,dx=1e-3,fval=None):
        if fval==None:
            fval = self.__call__(x)
        self._nGrad += 1
        f = self.__call__(x,True)
        df = self.__call__(x+dx,False)
        return fval,(df-f)/dx
    
    def get_taylor(self,x,dx=1e-3,fval=None):
        fval,grad = self.get_gradient(x,dx,fval)
        return Taylor1D(x,fval,grad)


class FunctionND:
    def __init__(self,func):
        self.func = func
        self._histXfull = 0
        self._histFfull = 0
        self._histXpart = 0
        self._histFpart = 0
        self._nEval = 0
        self._nGrad = 0
    
    def __call__(self,x,save=True):
        if self._nEval==0:
            self._nEval += 1
            fval = self.func(x)
            self._histXfull = np.array([x])
            self._histFfull = np.array([fval])
            if save:
                self._histXpart = np.array([x])
                self._histFpart = np.array([fval])
        else:
            for i,xval in enumerate(self._histXfull):
                if all(x==xval):
                    return self._histFfull[i]
            else:
                self._nEval += 1
                fval = self.func(x)
                self._histXfull = np.vstack([self._histXfull,x])
                self._histFfull = np.hstack([self._histFfull,fval])
                if save:
                    self._histXpart = np.vstack([self._histXpart,x])
                    self._histFpart = np.hstack([self._histFpart,fval])
        return fval

    def get_gradient(self,x,dx=1e-3,fval=None):
        if fval==None:
            fval = self.__call__(x,True)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i]+dx
            grad[i] = (self.__call__(X,False)-fval)/dx
        return fval, grad

    def get_taylor(self,x,dx=1e-3,fval=None):
        fval, grad = self.get_gradient(x,dx,fval)
        return TaylorND(x,fval,grad)


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
    func = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
    f = Function1D(func)
    approxF = f.get_taylor(0.4)
    
    plt.figure(1)
    plt.grid(True)
    plt.plot(x,[f(_x) for _x in x],'b-')
    plt.hold(True)
    plt.plot(x,[approxF(_x) for _x in x],'r-')
    plt.plot(f._histXpart,f._histFpart,'go')
    plt.show()

def test2d():
    func1d = lambda x: (5.0*x-2.0)**2.0*np.sin(12.*x-4.)
    func = lambda x: func1d(np.linalg.norm(x))
    f = FunctionND(func)
    approx = f.get_taylor(np.array([0.5,0.5]))
    print approx(np.array([0.5,0.5]))
    print approx(np.array([0.51,0.51]))
    print f._histFfull
    print f._histFpart

if __name__=="__main__":
    test2d()