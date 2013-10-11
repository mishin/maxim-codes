# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 19:35:10 2013

@author: Maxim
"""
from numpy import array, zeros, dot, copy, linspace, meshgrid, ravel
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class TaylorSeries1:
    def __init__(self,x0,f0,grad):
        self.x0 = array(x0)
        self.f0 = float(f0)
        self.grad = grad
    def __call__(self,x):
        return self.f0 + dot((array(x)-self.x0),self.grad)

class TestFunction():
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

    def get_gradient(self,x,dx=1e-10,fval=None):
        self.nGrad += 1
        if fval==None: fval = self(x)
        grad = zeros(len(x))
        for i in range(len(x)):
            X = copy(x)
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


class ScaledFunction():
    def __init__(self,funcLo,funcHi):
        self.funcHi = TestFunction(funcHi)
        self.funcLo = TestFunction(funcLo)
        self.nEval = 0
        self.xPrev = list()
        self.fPrev = list()
        self.betaPrev = list()
        self.taylor = True


def run_test_grad():
    def _fhigh(x):
        return 1.0/x[0] + 1.0/x[1] - 2.0
    
    x0 = array([3.0,3.0])
    func = TestFunction(_fhigh,'test function')
    ts = func.get_taylor(x0)
    
    x = linspace(0.5,10,30)
    y = linspace(0.5,10,30)
    X,Y = meshgrid(x,y)
    zs1 = array([func([xx,yy]) for xx,yy in zip(ravel(X), ravel(Y))])
    zs2 = array([ts([xx,yy]) for xx,yy in zip(ravel(X), ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    Z2 = zs2.reshape(X.shape)
    
    fig = plt.figure(1)
    ax = Axes3D(fig)
    ax.hold(True)
    ax.plot_wireframe(X,Y,Z1,color='y')
    ax.plot_wireframe(X,Y,Z2,color='b')
    ax.plot([x0[0]],[x0[1]],[func(x0)],'ro')
    plt.show()

if __name__=="__main__":
    run_test_grad()