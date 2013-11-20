# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 20:41:29 2013

@author: Maxim
"""

#from numpy import array, zeros, dot, copy, linspace, meshgrid, ravel, sin, pi
#from numpy.linalg import norm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import Rbf
from scipy.optimize import minimize
from miscTools import normalize, denormalize

class TaylorSeries1:
    def __init__(self,x0,f0,grad):
        self.x0 = np.array(x0)
        self.f0 = float(f0)
        self.grad = grad
    def __call__(self,x):
        return self.f0 + np.dot((np.array(x)-self.x0),self.grad)


class Function():
    def __init__(self,func,name=None):
        self.func = func
        self.name = name
        self.xHist = list()
        self.fHist = list()
        self.nEval = 0
        self.nGrad = 0

    def __call__(self,x,save=True):
        f = self.func(x)
        self.nEval += 1
        if save:
            self.xHist.append(x)
            self.fHist.append(f)
        return f

    def get_gradient(self,x,dx=1e-3,fval=None):
        self.nGrad += 1
        if fval==None: fval = self(x)
        if not hasattr(x,'__iter__'):
            grad = (self(x+dx) - fval)/dx
        else:
            grad = zeros(len(x))
            for i in range(len(x)):
                X = np.copy(x)
                X[i] = X[i]+dx
                grad[i] = (self(X,False)-fval)/dx
        return fval, grad

    def get_taylor(self,x,dx=1e-10,fval=None):
        fval, grad = self.get_gradient(x,dx,fval)
        return TaylorSeries1(x,fval,grad)

    def display(self,showHist=False):
        if self.name!=None:
            print self.name + '\n','-'*len(self.name)
        print 'Function evaluations: %d'%self.nEval
        print 'Gradient evaluations: %d'%self.nGrad
        if showHist:
            for xval,fval in zip(self.xHist,self.fHist):
                print xval, '| ',fval


class RbfMod():
    def __init__(self,x,y):
        if x.ndim>1:
            x = np.transpose(x)
            x = self._get_tuple(x)
            args = x + (y,)
        else:
            args = (x,y)
        print args
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


def run_test1():
    x = np.array([[0,0],[0,1],[1,0],[1,1],[0.5,0.5]])
    func = lambda x: x[0]**2 + 0.5*x[1]
    
    y = np.array([func(_x) for _x in x])
    approx = RbfMod(x,y)
    dx = 0.05
    xmsh = ymsh = np.arange(0,1,dx)
    Xmsh, Ymsh = np.meshgrid(xmsh,ymsh)
    Zmsh = np.array([approx(np.array([_x,_y])) for _x,_y in zip(np.ravel(Xmsh),np.ravel(Ymsh))])
    Zmsh = Zmsh.reshape(Xmsh.shape)
    Zmsh2 = np.array([func(np.array([_x,_y])) for _x,_y in zip(np.ravel(Xmsh),np.ravel(Ymsh))])
    Zmsh2 = Zmsh2.reshape(Xmsh.shape)
    fig = plt.figure(1)
    ax = Axes3D(fig)
    ax.hold(True)
    ax.plot_wireframe(Xmsh,Ymsh,Zmsh)
    ax.plot_wireframe(Xmsh,Ymsh,Zmsh2,color='r')
    plt.show()

def run_test2():
    x = np.array([0,0.1,0.5,0.8,1.0])
    func = lambda x: (5.0*x-2.0)**2.0*np.sin(12.*x-4.)
    y = np.array([func(_x) for _x in x])
    approx = RbfMod(x,y)
    xx = np.linspace(0,1)
    y2 = np.array([approx(_x) for _x in xx])
    y3 = np.array([func(_x) for _x in xx])

    plt.show()
    plt.hold(True)
    plt.plot(xx,y2,'b-')
    plt.plot(xx,y3,'r-')
    plt.plot(x,y,'ko')
    plt.show()
    
if __name__=="__main__":
    #run_test1()
    run_test2()