# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 19:35:10 2013

@author: Maxim
"""
from numpy import array, zeros, dot, copy, linspace, meshgrid, ravel
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import Rbf
from scipy.optimize import minimize

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
        if not hasattr(x,'__iter__'):
            grad = (self(x+dx) - fval)/dx
        else:
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
    def __init__(self,funcLo,funcHi,minXnum=3):
        self.funcHi = TestFunction(funcHi)
        self.funcLo = TestFunction(funcLo)
        self.nEval = 0
        self.xPrev = list()
        self.fPrev = list()
        self.betaPrev = list()
        self.nMin = minXnum
        self.oneDim = False

    def construct_scaling_model(self,x0,f0=None):
        if hasattr(x0,'__iter__'):
            self.oneDim = True
        self.x0 = x0
        if f0==None:
            f0 = self.funcHi(x0)
        self.xPrev.append(x0)
        self.fPrev.append(f0)
        if len(self.xPrev)<=self.nMin:
            fHi, gradHi = self.funcHi.get_gradient(x0,fval=f0)
            fLo, gradLo = self.funcLo.get_gradient(x0)
            betaGrad = (gradHi*fLo - gradLo*fHi) / (fLo**2)
            beta = fHi/fLo
            self.betaPrev.append(beta)
            self.beta = TaylorSeries1(x0,beta,betaGrad)
        else:
            fLo = self.funcLo(x0)
            beta = f0/fLo
            self.betaPrev.append(beta)
            if self.oneDim:
                xPrev = array(self.xPrev)
            else:
                m = len(self.xPrev[0])
                xPrev = tuple()
                xnew = zeros(m)
                for i in range(m):
                    for xval in self.xPrev:
                        xnew[i] = xval[i]
                    xPrev = xPrev + (xnew,)
            xrbf = tuple()
            xrbf = xrbf + (xPrev,self.betaPrev)
            self.beta = Rbf(*xrbf)

    def __call__(self,x):
        self.nEval += 1
        return self.beta(x) * self.funcLo(x)
    
    def get_thrust_region_ratio(self,x):
        fHi = self.funcHi(x)
        fScaled = self(x)
        fScaled0 = self(self.x0)
        if fScaled==fScaled0:
            return float('inf')
        else:
            rho = (fScaled0 - fHi) / (fScaled0 - fScaled)
            return rho

