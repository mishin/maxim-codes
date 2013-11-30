# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 01:20:13 2013

@author: Maxim
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class FunctionND:
    """ Class handling multidimensional functions and provides additional 
    functionality to reduce computational time for costly functions and provides
    easy tools for evaluation number of iterations, and history of calls
    """
    def __init__(self,func):
        """
        Parameters
        ----------
        
        func: function
            callable function
        """
        self.func = func
        self._histXfull = 0
        self._histFfull = 0
        self._histXpart = 0
        self._histFpart = 0
        self._histIsPart = list()
        self._histTolerance = 1e-15
        self._nEval = 0
        self._nGrad = 0
    
    def __call__(self,x,save=True):
        if not hasattr(x,'__iter__'):
            x = np.array([x])
        if self._nEval==0:
            self._nEval += 1
            fval = self._feval(x)
            self._histXfull = np.array([x])
            self._histFfull = np.array([fval])
            if save:
                self._histIsPart.append(1)
                self._histXpart = np.array([x])
                self._histFpart = np.array([fval])
            else:
                self._histIsPart.append(0)
        else:
            for i,xval in enumerate(self._histXfull):
                if all(abs(x-xval)<=self._histTolerance):
                    if save:
                        self._save_part_history(self._histXfull[i],self._histFfull[i])
                    return self._histFfull[i]
            else:
                self._nEval += 1
                fval = self._feval(x)
                self._histXfull = np.vstack([self._histXfull,x])
                self._histFfull = np.hstack([self._histFfull,fval])
                if save:
                    self._histIsPart.append(1)
                    self._histXpart = np.vstack([self._histXpart,x])
                    self._histFpart = np.hstack([self._histFpart,fval])
                else:
                    self._histIsPart.append(0)
        return fval
    
    def _save_part_history(self,x,fval):
        exist = False
        for xval in self._histXpart:
            if xval==x:
                exist = True
                break
        if not exist:
            self._histIsPart.append(1)
            self._histXpart = np.vstack([self._histXpart,x])
            self._histFpart = np.hstack([self._histFpart,fval])

    def _feval(self,x):
        if len(x)==1:
            return self.func(x[0])
        else:
            return self.func(x)

    def get_gradient(self,x,dx=1e-3,fval=None):
        """
        Returns gradient at point x evaluated using forward difference scheme.
        
        Parameters
        ----------
        
        x: 1d array
            point at which gradient will be evaluated
        dx: float
            delta x for numerical calculation of gradient
        fval : float
            function value at given point. If fval==None then function will be 
            evaluated.
        """
        if not hasattr(x,'__iter__'):
            x = np.array([x])
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
    
    def write_history(self,filePath):
        fid = open(filePath,'wt')
        for isPart,x,f in zip(self._histIsPart,self._histXfull,self._histFfull):
            fid.write('%d\t'%isPart)
            for xval in x:
                fid.write('%.15f\t'%xval)
            fid.write('%.15f\n'%f)
        fid.close()
    
    def read_history(self,filePath):
        fid = open(filePath,'rt')
        lines = fid.readlines()
        fid.close()
        self._initialize_history()
        for i,line in enumerate(lines):
            if line.strip()!='':
                segLine = line.split()
                isPart = int(segLine[0])
                fval = float(segLine[-1])
                xval = np.array([float(val) for val in segLine[1:-1]])
                self._histIsPart.append(isPart)
                self._histFfull = np.hstack([self._histFfull,fval])
                if i==0:
                    self._histXfull = xval
                    if bool(isPart):
                        self._histFpart = np.hstack([self._histFpart,fval])
                        self._histXpart = xval
                else:
                    self._histXfull = np.vstack([self._histXfull,xval])
                    if bool(isPart):
                        self._histFpart = np.hstack([self._histFpart,fval])
                        if self._histXpart==list():
                            self._histXpart = xval
                        else:
                            self._histXpart = np.vstack([self._histXpart,xval])

    def _initialize_history(self):
        self._histXfull = list()
        self._histFfull = list()
        self._histXpart = list()
        self._histFpart = list()
        self._histIsPart = list()


class TaylorND:
    def __init__(self,x,f,gradF):
        self.x0 = x
        self.f0 = f
        self.grad = gradF
    def __call__(self,x):
        return self.f0 + np.dot((np.array(x)-self.x0),self.grad)


def forrester(x):
    return (6.0*x-2.0)**2.0*np.sin(12.*x-4.)

def test_1():
    x = np.linspace(0,1,50)
    func = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
    f = FunctionND(func)
    approxF = f.get_taylor(0.4,1e-10)
    
    plt.figure(1)
    plt.grid(True)
    plt.plot(x,[f(_x) for _x in x],'b-')
    plt.hold(True)
    plt.plot(x,[approxF(_x) for _x in x],'r-')
    plt.plot(f._histXpart,f._histFpart,'go')
    f.write_history('tmp_history77.txt')
    f.read_history('tmp_history77.txt')
    plt.show()

def test_2():
    func = lambda x: x[0]+x[1]-x[1]**2
    f = FunctionND(func)
    print f(np.array([0,4]))
    print f(np.array([5,4]),False)
    print f._nEval
    print f._histFfull
    print f._histFpart

if __name__=="__main__":
    test_2()