# -*- coding: utf-8 -*-
"""
Created on Wed Dec 04 16:06:11 2013

@author: Maxim
"""

import numpy as np


class Function:
    def __init__(self,func):
        self.func        = func
        self._histXfull  = list()
        self._histXpart  = list()
        self._histFfull  = list()
        self._histFpart  = list()
        self._histIsPart = list()
        self._histTol    = 1e-15
        self._dx         = 1e-6
        self._nEval      = 0
        self._nGrad      = 0
    
    def __call__(self,x,save=True):
        if not hasattr(x,'__iter__'):
            x = np.array([x])
        if self._nEval==0:
            fval = self._call_first_iteration(x,save)
        else:
            fval = self._call_next_iteration(x,save)
        return fval
    
    def _call_first_iteration(self,x,save=True):
        self._nEval += 1
        fval = self._feval(x)
        self._histXfull = np.array([x])
        self._histFfull = np.array([fval])
        if save:
            self._histIsPart.append(1)
            self._histXpart = self._histXfull[-1]
            self._histFpart = self._histFfull[-1]
        else:
            self._histIsPart.append(0)
        return fval
    
    def _call_next_iteration(self,x,save=True):
        for i,xval in enumerate(self._histXfull):
            if all(abs(x-xval)<=self._histTol):
                return self._histFfull[i]
        else:
            self._nEval += 1
            fval = self._feval(x)
            self._histXfull = np.vstack([self._histXfull, x])
            self._histFfull = np.hstack([self._histFfull, fval])
            if save:
                self._histIsPart.append(1)
                self._histXpart = np.vstack([self._histXpart,self._histXfull[-1]])
                self._histFpart = np.hstack([self._histFpart,self._histFfull[-1]])
            else:
                self._histIsPart.append(0)
            return fval

    def _feval(self,x):
        if len(x)==1:
            return self.func(x[0])
        else:
            return self.func(x)
    
    def get_gradient(self,x,dx=None):
        if not hasattr(x,'__iter__'):
            x = np.array([x])
        if dx==None:
            dx = self._dx
        fval = self.__call__(x,True)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i]+dx
            grad[i] = (self.__call__(X,False)-fval)/dx
        return fval,grad

    def write_history(self,filePath):
        fid = open(filePath,'wt')
        for isPart,x,f in zip(self._histIsPart,self._histXfull,self._histFfull):
            fid.write('%d\t'%isPart)
            for xval in x:
                fid.write('%.15f\t'%xval)
            fid.write('%.15f\n'%f)
        fid.close()

    def get_taylor(self,x,dx=1e-3,fval=None):
        fval, grad = self.get_gradient(x,dx)
        return Taylor1(x,fval,grad)

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


class Taylor1:
    def __init__(self,x,f,gradF):
        self.x0 = x
        self.f0 = f
        self.grad = gradF
    def __call__(self,x):
        return self.f0 + np.dot((np.array(x)-self.x0),self.grad)


# --- Debug section ---
def forrester(x):
    return (6.0*x-2.0)**2.0*np.sin(12.*x-4.)

def test_1():
    x = np.linspace(0,1,50)
    func = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
    func2 = lambda x: (5*x-2.)**2*np.sin(12*x-4)
    f = Function(func)
    approxF = f.get_taylor(0.4)
    
    plt.figure(1)
    plt.grid(True)
    plt.plot(x,[f(_x) for _x in x],'b-')
    plt.hold(True)
    plt.plot(x,[approxF(_x) for _x in x],'r-')
    plt.plot(f._histXpart,f._histFpart,'ro')
    plt.plot(x,func2(x),'k-')
    f.write_history('tmp_history77.txt')
    f.read_history('tmp_history77.txt')
    plt.show()

def test_2():
    func = lambda x: x[0]+x[1]-x[1]**2
    f = Function(func)
    print f(np.array([0,4]))
    print f(np.array([5,4]),False)
    print f._nEval
    print f._histFfull
    print f._histFpart

if __name__=="__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    test_1()