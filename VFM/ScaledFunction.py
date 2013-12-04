# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 14:02:57 2013

@author: Maxim
"""
from functionHandling import *
from scipy.interpolate import Rbf

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


class HybridScaling:
    """
    weight=0 -> additive
    """
    def __init__(self,fHigh,fLow,weight=0.5,warmUpRun=3,dx=1e-4):
        self.w = weight
        self.fAdd = ScaledFunction(fHigh,fLow,'add',warmUpRun,dx)
        self.fMult = ScaledFunction(fHigh,fLow,'mult',warmUpRun,dx)
    
    def construct_scaling_model(self,x0):
        self.x0 = x0
        if self.w==0:
            self.fAdd.construct_scaling_model(x0)
            self.f0 = self.fAdd.fH0
        elif self.w==1:
            self.fMult.construct_scaling_model(x0)
            self.f0 = self.fMult.fH0
        else:
            self.fAdd.construct_scaling_model(x0)
            self.f0 = self.fAdd.fH0
            self.fMult.construct_scaling_model(x0)
    
    def initialize_by_points(self,X):
        if self.w==0:
            self.fAdd.initialize_by_points(X)
        elif self.w==1:
            self.fMult.initialize_by_points(X)
        else:
            self.fAdd.initialize_by_points(X)
            self.fMult.initialize_by_points(X)
    
    def get_trust_region_ratio(self,x):
        fSc = self.__call__(x)
        if self.w==0:
            fH = self.fAdd.fHigh(x,True)
            fL = self.fAdd.fLow(x,True)
        else:
            fH = self.fMult.fHigh(x,True)
            fL = self.fMult.fLow(x,True)
        fSc0 = self.f0
        if fSc==fSc0:
            return float('inf')
        else:
            return (fSc0 - fH) / (fSc0 - fSc)
    
    def __call__(self,x,save=False):
        if self.w==0:
            return self.fAdd(x,save)
        elif self.w==1:
            return self.fMult(x,save)
        else:
            return self.w*self.fMult(x,save) + (1.-self.w)*self.fAdd(x,save)

class ScaledFunction:
    """
    fhigh, flow - functions
    fH, fL - function values
    """
    def __init__(self,fHigh,fLow,scaling='add',warmUpRun=3,dx=1e-4):
        self.fHigh = FunctionND(fHigh)
        self.fLow = FunctionND(fLow)
        self.scalingType = scaling
        self.warmUpRun = warmUpRun
        self._dx = dx
        self._histScaling = list()

    def construct_scaling_model(self,x0):
        fH = self.fHigh(x0,True)
        fL = self.fLow(x0,True)
        self.x0 = x0
        self.fH0 = fH
        scalingFactor = self._get_scaling_factor(fH,fL)
        _list = zip(self.fHigh._histFpart,self.fLow._histFpart)
        _hist = [self._get_scaling_factor(_fH,_fL) for _fH,_fL in _list]
        self._histScaling = np.array([_hist])
        if len(self.fHigh._histFpart)<self.warmUpRun:
            scalingGrad = self._get_scaling_factor_gradient(x0,fH,fL)
            self.scalingFunc = TaylorND(x0,scalingFactor,scalingGrad)
        else:
            self.scalingFunc = RbfMod(self.fHigh._histXpart,self._histScaling)

    def __call__(self,x,save=False):
        if self.scalingType=='add':
            return self.scalingFunc(x) + self.fLow(x,save)
        elif self.scalingType=='mult':
            return self.scalingFunc(x) * self.fLow(x,save)
    
    def _get_scaling_factor(self,fH,fL):
        if self.scalingType=='add':
            return fH-fL
        elif self.scalingType=='mult':
            return fH/fL
        
    def _get_scaling_factor_gradient(self,x,fH,fL):
        fH, gradH = self.fHigh.get_gradient(x,self._dx,fH)
        fL, gradL = self.fLow.get_gradient(x,self._dx,fL)
        if self.scalingType=='add':
            return gradH - gradL
        elif self.scalingType=='mult':
            return (gradH*fL - gradL*fH) / (fL*fL)

    def get_trust_region_ratio(self,x):
        fSc = self.__call__(x)
        fH = self.fHigh(x,True)
        fL = self.fLow(x,True)
        fSc0 = self.fH0
        if fSc==fSc0:
            return float('inf')
        else:
            return (fSc0 - fH) / (fSc0 - fSc)
    
    def initialize_by_points(self,X):
        for x in X[:-1]:
            self.fHigh(x,True)
            self.fLow(x,True)
        self.construct_scaling_model(X[-1])
    
    def initialize_from_file(self,fileHigh,fileLow=None):
        pass
    
    def get_gradient(self,x):
        fval = self.__call__(x,False)
        grad = np.zeros(len(x))
        for i in range(len(x)):
            X = np.copy(x)
            X[i] = X[i] + self._dx
            grad[i] = (self.__call__(X,False)-fval)/self._dx


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
        elif self.eta1<rho<self.eta2:
            delta = delta0
        else:
            delta = delta0 * self._get_gamma(x0,xnew,delta0)
        if rho>0:
            x0new = xnew
        else:
            x0new = x0
        self.delta0 = delta
        return delta, x0new
    
    def _get_gamma(self,x0,xnew,delta):
        err = np.linalg.norm(x0-xnew)
        if hasattr(x0,'__iter__'):
            err0 = np.linalg.norm(np.ones(len(x0))*delta)
        else:
            err0 = delta
        if err==err0:
            return self.c2
        else:
            return 1.0



def run_test1():
    x = np.linspace(0,1,50)
    fhigh = forrester
    flow = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5.
    
    fsc = ScaledFunction(fhigh,flow,'add',2,True)
    fsc.construct_scaling_model(0.1)
    fsc.construct_scaling_model(0.9)
    fsc.construct_scaling_model(0.5)
    print fsc.get_trust_region_ratio(0.8)
    print fsc.fHigh._nEval
    fsc.construct_scaling_model(0.8)
    print fsc.fHigh._nEval

    y1 = np.array([fhigh(_x) for _x in x])
    y2 = np.array([flow(_x) for _x in x])
    y3 = np.array([fsc(_x) for _x in x])
    plt.figure(1)
    plt.hold(True)
    plt.plot(x,y1,'r-')
    plt.plot(x,y2,'b-')
    plt.plot(x,y3,'g-')
    plt.plot(fsc.fHigh._histXpart, fsc.fHigh._histFpart,'ko')
    plt.show()

def run_test2():
    fhigh = lambda x: forrester(np.linalg.norm(x))
    flow = lambda x: fhigh(x)+20*np.linalg.norm(x) + 50*sum(x)
    
    fsc = HybridScaling(fhigh,flow,0.7,3)
    fsc.construct_scaling_model(np.array([0.9,0.9]))
    fsc.construct_scaling_model(np.array([0.5,0.5]))
    fsc.construct_scaling_model(np.array([0.4,0.8]))
    fsc.construct_scaling_model(np.array([0.1,0.1]))

    lb = [0, 0.0]
    ub = [1., 1.]
    dx = 0.05
    x = y = np.arange(lb[0],ub[0],dx)
    X, Y = np.meshgrid(x,y)
    Z1 = np.array([fsc.fMult.fHigh(np.array([_x,_y]),False) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z2 = np.array([fsc.fMult.fLow(np.array([_x,_y]),False) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z3 = np.array([fsc(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    
    Z1 = Z1.reshape(X.shape)
    Z2 = Z2.reshape(X.shape)
    Z3 = Z3.reshape(X.shape)
    
    fig = plt.figure(1)
    ax1 = Axes3D(fig)
    ax1.hold(True)
    ax1.plot_wireframe(X,Y,Z1,color='b')
    ax1.plot_wireframe(X,Y,Z2,color='r')
    ax1.plot_wireframe(X,Y,Z3,color='y')
    ax1.plot(fsc.fMult.fHigh._histXpart[:,0],fsc.fMult.fHigh._histXpart[:,1],fsc.fMult.fHigh._histFpart,'ro')
    plt.show()

if __name__=="__main__":
    run_test2()