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
        self.fHist = np.array([])
        self.nEval = 0
        self.nGrad = 0

    def __call__(self,x,save=True):
        if hasattr(x,'__iter__'):
            self.nvar = len(x)
        else:
            self.nvar = 1
        
        if self.nEval==0:
            idx=None
        else:
            idx = self._check(x)
            
        if idx==None:
            f = self.func(x)
            if save:
                self.nEval += 1
                self._save_history(x,f)
        else:
            f = self.fHist[idx]
        return f
    
    def _save_history(self,x,f):
        if self.nEval==1:
            self.xHist = np.array([[x]])
        else:
            self.xHist = np.vstack([self.xHist, x])
        self.fHist = np.append(self.fHist,f)

    def _check(self,x):
        for i,xh in enumerate(self.xHist):
            if all(x==xh):
                return i
                break
        else:
            return None

    def get_gradient(self,x,dx=1e-3,fval=None):
        self.nGrad += 1
        if fval==None:
            fval = self(x)
        if not hasattr(x,'__iter__'):
            grad = (self(x+dx,False) - fval)/dx
        else:
            grad = np.zeros(len(x))
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


class ScaledFunction:
    def __init__(self,funcLow,funcHigh,warmUpRun,scalingType='add'):
        self.funcLow  = Function(funcLow)
        self.funcHigh = Function(funcHigh)
        self.nEval = 0
        self.historyScaledFunc = list()
        self.historyScaledFact = list()
        self.nWarmUp = warmUpRun
        self.type = scalingType
        self._dx = 1e-3

    def construct_scaling_model(self,x0,fHighNew=None,fLowNew=None):
        if fHighNew==None:
            fHighNew = self.funcHigh(x0)
        if fLowNew==None:
            fLowNew = self.funcLow(x0)
        if self.funcHigh.nEval<=self.nWarmUp:
            print 'taylor'
            self._construct_taylor_series(x0,fHighNew,fLowNew)
        else:
            print 'rbf'
            self._construct_rbf(x0,fHighNew,fLowNew)
    
    def _construct_taylor_series(self,x,fHigh,fLow):
        fHigh, gradHigh = self.funcHigh.get_gradient(x,fval=fHigh)
        fLow,  gradLow  = self.funcLow.get_gradient(x,fval=fLow)
        if self.type=='add':
            gradScaling0 = (gradHigh - gradLow)
        elif self.type=='mult':
            gradScaling0 = (gradHigh*fLow - gradLow*fHigh)/(fLow**2)
        scFactor0 = self._get_scaling_factor(x,fHigh,fLow)
        self.scalingFactor = TaylorSeries1(x,scFactor0,gradScaling0)
        self.historyScaledFact.append(scFactor0)
    
    def _construct_rbf(self,x,fHigh,fLow):
        scFactor0 = self._get_scaling_factor(x,fHigh,fLow)
        self.historyScaledFact.append(scFactor0)
        self.scalingFactor = RbfMod(self.funcHigh.xHist,self.historyScaledFact)

    def _get_scaling_factor(self,x,fHigh,fLow):
        if self.type=='add':
            return fHigh-fLow
        elif self.type=='mult':
            return fHigh/fLow
    
    def __call__(self,x):
        self.nEval += 1
        fLow = self.funcLow(x)
        if self.type=='add':
            return self.scalingFactor(x) + fLow
        elif self.type=='mult':
            return fLow*self.scalingFactor(x)
    
    def get_trust_region_ratio(self,x,fHigh=None):
        fHigh0 = self.funcHigh.fHist[-1]
        if fHigh==None:
            fHigh = self.funcHigh(x)
        fScaled = self.__call__(x)
        fscaled0 = fHigh0
        if fScaled==fscaled0:
            return float('inf'),fHigh
        else:
            rho = (fHigh0 - fHigh) / (fscaled0 - fScaled)
            return rho, fHigh
    
        
class TrustRegionManagement:
    def __init__(self,delta,eta1=0.25,eta2=0.75,eta3=1.25,c1=0.3,c2=2.0):
        self.eta1 = eta1
        self.eta2 = eta2
        self.eta3 = eta3
        self.c1 = c1
        self.c2 = c2
        self.delta0 = delta
    
    def adjust(self,rho,x0,xnew):
        delta0 = self.delta0
        if rho<=self.eta1 or rho>=self.eta3:
            delta = self.c1*delta0
            x0new = x0
        elif self.eta1<rho<self.eta2:
            delta = delta0
            x0new = xnew
        else:
            x0new = xnew
            xerr = np.array(xnew-x0)
            if xerr.any==delta0:
                delta = self.c2*delta0
            else:
                delta = delta0
        self.delta0 = delta
        return delta, x0new

def run_test3():
    func1 = lambda x: (1.-x[0])**2 + 100*(x[1]-x[0]**2)**2
    func = Function(func1)
    
    func(np.array([0,0]))
    func(np.array([5,1]))
    func(np.array([2,3]))
    
    rbf = RbfMod(func.xHist,func.fHist)

def run_test4():
    func1 = lambda x: (6.*x-2.)**2.*np.sin(12.*x-4.)
    func = Function(func1)
    func(0)
    func(0.5)
    func(0.0)
    func(-1)

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
    #x = np.array([0,0.1,0.5,0.8,1.0])
    x = np.linspace(0,1,8)
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
    run_test4()