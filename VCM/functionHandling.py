# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 01:20:13 2013

@author: Maxim
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Function1D:
    def __init__(self,func):
        self.func = func
        self._histXfull = list()
        self._histFfull = list()
        self._histXpart = list()
        self._histFpart = list()
        self._histIsPart = list()
        self._nEval = 0
        self._nGrad = 0
        self._histTolerance = 1e-15
    
    def __call__(self,x,save=True):
        for i,xval in enumerate(self._histXfull):
            if abs(x-xval)<=self._histTolerance:
                return self._histFfull[i]
        else:
            self._nEval += 1
            f = self.func(x)
            self._histXfull = np.hstack([self._histXfull,x])
            self._histFfull = np.hstack([self._histFfull,f])
            if save:
                self._histIsPart.append(1)
                self._histXpart = np.hstack([self._histXpart,x])
                self._histFpart = np.hstack([self._histFpart,f])
            else:
                self._histIsPart.append(0)
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
    
    def write_history(self,filePath):
        """ writes full and part history to a file specified for later use
        """
        fid = open(filePath,'wt')
        for isPart,x,f in zip(self._histIsPart,self._histXfull,self._histFfull):
            fid.write('%d\t%.15f\t%.15f\n'%(isPart,x,f))
        fid.close()
    
    def read_history(self,filePath):
        """ reads in history file in format of write_history method.
        Note: previous history is deleted
        """
        fid = open(filePath,'rt')
        lines = fid.readlines()
        fid.close()
        self._initialize_history()
        for i,line in enumerate(lines):
            if line.strip()!='':
                segLine = line.split()
                isPart = int(segLine[0])
                xval = float(segLine[1])
                fval = float(segLine[2])
                self._histXfull.append(xval)
                self._histFfull.append(fval)
                self._histIsPart.append(isPart)
                if bool(isPart):
                    self._histXpart.append(xval)
                    self._histFpart.append(fval)
        self._histFfull = np.array(self._histFfull)
        self._histXfull = np.array(self._histXfull)
        self._histFpart = np.array(self._histFpart)
        self._histXpart = np.array(self._histXpart)
    
    def _initialize_history(self):
        self._histXfull = list()
        self._histFfull = list()
        self._histXpart = list()
        self._histFpart = list()
        self._histIsPart = list()


class FunctionND:
    def __init__(self,func):
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
        if self._nEval==0:
            self._nEval += 1
            fval = self.func(x)
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
                    return self._histFfull[i]
            else:
                self._nEval += 1
                fval = self.func(x)
                self._histXfull = np.vstack([self._histXfull,x])
                self._histFfull = np.hstack([self._histFfull,fval])
                if save:
                    self._histIsPart.append(1)
                    self._histXpart = np.vstack([self._histXpart,x])
                    self._histFpart = np.hstack([self._histFpart,fval])
                else:
                    self._histIsPart.append(0)
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


def forrester(x):
    return (6.0*x-2.0)**2.0*np.sin(12.*x-4.)

def test1d():
    x = np.linspace(0,1,4)
    func = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
    f = Function1D(func)
    approxF = f.get_taylor(0.4)
    
    plt.figure(1)
    plt.grid(True)
    plt.plot(x,[f(_x) for _x in x],'b-')
    plt.hold(True)
    plt.plot(x,[approxF(_x) for _x in x],'r-')
    plt.plot(f._histXpart,f._histFpart,'go')
    f.write_history('tmp_history77.txt')
    f.read_history('tmp_history77.txt')
    plt.show()

def test2d():
    func1d = lambda x: (5.0*x-2.0)**2.0*np.sin(12.*x-4.)
    func = lambda x: func1d(np.linalg.norm(x))
    f = FunctionND(func)
    approx = f.get_taylor(np.array([0.7,0.7]),1e-10)
    lb = [0, 0.0]
    ub = [1., 1.]
    dx = 0.05
    x = y = np.arange(lb[0],ub[0],dx)
    X, Y = np.meshgrid(x,y)
    Z = np.array([f(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z = Z.reshape(X.shape)
    Z1 = np.array([approx(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z1 = Z1.reshape(X.shape)
    
    fig = plt.figure(1)
    ax1 = Axes3D(fig)
    ax1.plot_wireframe(X,Y,Z)
    ax1.hold(True)
    ax1.plot_wireframe(X,Y,Z1,color='r')
    f.write_history('tmp_history88.txt')
    f.read_history('tmp_history88.txt')
    plt.show()

def test2dsave():
    func1d = lambda x: (5.0*x-2.0)**2.0*np.sin(12.*x-4.)
    func = lambda x: func1d(np.linalg.norm(x))
    
    f = FunctionND(func)
    f(np.array([0,5]))
    f(np.array([1,5]))
    f(np.array([3,4]))
    f(np.array([9,6]))
    f.get_gradient(np.array([0.5,0.6]))
    f.write_history('tmp_history88.txt')
    print f._histFpart
    print f._histFfull
    f.read_history('tmp_history88.txt')
    print f._histFpart
    print f._histFfull
if __name__=="__main__":
    test2d()