# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 19:35:10 2013

@author: Maxim
"""
from numpy import array, zeros, dot, copy, linspace, meshgrid, ravel, sin, pi
from numpy.linalg import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import Rbf
from scipy.optimize import minimize
from miscTools import normalize, denormalize

class TrustRegionManagement():
    def __init__(self,delta,eta1=0.25,eta2=0.75,eta3=1.25,c1=0.3,c2=2.0):
        self.eta1 = eta1
        self.eta2 = eta2
        self.eta3 = eta3
        self.c1 = c1
        self.c2 = c2
        self.deltaOld = delta
    def adjust(self,rho,err):
        if rho<=self.eta1 or rho>=self.eta3:
            deltaNew = self.deltaOld *self.c1
        elif self.eta2< rho <self.eta3:
            deltaNew = self.deltaOld *self.c2
#            if err==self.deltaOld:
#                #FIXME: there is no need in keeping same delta if approximation is good
#                deltaNew = self.deltaOld *self.c2
#            else:
#                deltaNew = self.deltaOld
        else:
            deltaNew = self.deltaOld
        self.deltaOld = deltaNew
        return deltaNew

class RbfMod():
    def __init__(self,args):
        self.rbf = Rbf(*args)
    def __call__(self,x):
        if hasattr(x,'__iter__'):
            _x = tuple()
            for xx in x:
                _x += (xx,)
            return self.rbf(*_x)
        else:
            return self.rbf(x)


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

    def get_gradient(self,x,dx=0.01,fval=None):
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


class HybridScaledFunction():
    """
    weight = 1.0 - multiplicative
    weight = 0.0 - additive
    """
    def __init__(self,funcLo,funcHi,minXnum=4,weight=0.5):
        self.w = weight
        self.nEval = 0
        self.xPrev = list()
        self.fPrev = list()
        self.fAdd  = ScaledFunction(funcLo,funcHi,minXnum,'add')
        self.fMult = ScaledFunction(funcLo,funcHi,minXnum,'mult')
        self.funcHi = TestFunction(funcHi)
        self.funcLo = TestFunction(funcLo)
    
    def __call__(self,x):
        self.nEval += 1
        return self.w*self.fMult(x) + (1.-self.w)*self.fAdd(x)
    
    def construct_scaling_model(self,x0,f0=None):
        self.x0 = x0
        if f0==None:
            f0 = self.fAdd.funcHi(x0)
        self.xPrev.append(x0)
        self.fPrev.append(f0)
        self.fAdd.construct_scaling_model(x0,f0)
        self.fMult.construct_scaling_model(x0,f0)
    
    def _initialize_by_doe_points(self,xnew,fnew=None):
        if fnew==None:
            fnew = list()
            for i,xx in enumerate(xnew):
                fnew.append(self.fAdd.funcHi(xx))
        for _f,_x in zip(fnew,xnew):
            self.fPrev.append(_f)
            self.xPrev.append(_x)
        self.fAdd._initialize_by_doe_points(xnew,fnew)
        self.fMult._initialize_by_doe_points(xnew,fnew)

    def get_trust_region_ratio(self,x,fHi=None):
        if fHi==None:
            fHi = self.fAdd.funcHi(x)
        fScaled = self(x)
        fScaled0 = self(self.x0)
        if fScaled==fScaled0:
            return float('inf'),fHi
        else:
            rho = (fScaled0 - fHi) / (fScaled0 - fScaled)
            return rho,fHi
    
    def display(self):
        self.funcHi.display()
        self.funcLo.display()

class ScaledFunction():
    def __init__(self,funcLo,funcHi,minXnum=4,scalingType='mult'):
        self.funcHi = TestFunction(funcHi)
        self.funcLo = TestFunction(funcLo)
        self.nEval = 0
        self.xPrev = list()
        self.fPrev = list()
        self.betaPrev = list()
        self.gammaPrev = list()
        self.nMin = minXnum
        self.oneDim = False
        self.type = scalingType
        self.dx = 1e-3

    def construct_scaling_model(self,x0,f0=None):
        if hasattr(x0,'__iter__'):
            self.oneDim = False
        else:
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
            #FIXME beta calculation here
            beta = self.get_beta(x0,fHi)
            self.betaPrev.append(beta)
            self.beta = TaylorSeries1(x0,beta,betaGrad)
        else:
            beta = self.get_beta(x0,f0)
            self.betaPrev.append(beta)
            xrbf = tuple()
            if self.oneDim:
                xrbf = xrbf + (array(self.xPrev),)
            else:
                n = len(self.xPrev); m = len(self.xPrev[0])
                for i in range(m):
                    xnew = zeros(n)
                    for j,xval in enumerate(self.xPrev):
                        xnew[j] = xval[i]
                    xrbf = xrbf + (xnew,)
            xrbf = xrbf + (self.betaPrev,)
            self.beta = RbfMod(xrbf)
    
    def _initialize_by_doe_points(self,xnew,fnew=None,flow=None):
        print '--> initializing points'
        if fnew==None:
            fnew = list()
            for i,xx in enumerate(xnew):
                fnew.append(self.funcHi(xx))
        if flow==None:
            for _f,_x in zip(fnew,xnew):
                self.fPrev.append(_f)
                self.xPrev.append(_x)
                beta = self.get_beta(_x,_f)
                self.betaPrev.append(beta)
        else:
            for _f,_x,_flow in zip(fnew,xnew,flow):
                self.fPrev.append(_f)
                self.xPrev.append(_x)
                beta = self.get_beta(_x,_f,_flow)
                self.betaPrev.append(beta)
        print '\n--> initialization completed'

    def __call__(self,x):
        self.nEval += 1
        if self.type=='add':
            return self.beta(x) + self.funcLo(x)
        elif self.type=='mult':
            return self.beta(x) * self.funcLo(x)
        
    def derivative(self,x,dx=None):
        if dx==None:
            dx = self.dx
        grad = zeros(len(x))
        fval = self.__call__(x)
        for i,xx in enumerate(x):
            X = copy(x)
            X[i] = X[i]+dx
            grad[i] = (self.__call__(X)-fval)/dx
        return grad
    
    def get_trust_region_ratio(self,x,fHi=None):
        if fHi==None:
            fHi = self.funcHi(x)
        fScaled = self(x)
        fScaled0 = self(self.x0)
        if fScaled==fScaled0:
            return float('inf'),fHi
        else:
            rho = (fScaled0 - fHi) / (fScaled0 - fScaled)
            return rho,fHi
    
    def get_beta(self,x,fHi=None,fLow=None):
        if fHi==None:
            fHi = self.funcHi(x)
        if fLow==None:
            fLow = self.funcLo(x)
        if self.type=='add':
            return fHi - fLow
        elif self.type=='mult':
            return fHi / fLow

class VCMoptimization:
    """variables will be normalized
    """
    def __init__(self,eta1=0.25,eta2=0.75,eta3=1.25,c1=0.25,c2=2.0,delta=1.0):
        self.eta1 = float(eta1)
        self.eta2 = float(eta2)
        self.eta3 = float(eta3)
        self.c1   = float(c1)
        self.c2   = float(c2)
        self.delta= float(delta)
        self.constr    = list()
        self.constrVcm = list()
        self.nVcmConstr = 0
        self.tolX = 1.0e-6
        self.tolF = 1.0e-6
        self.iterMax = 50
        self.vcmConstr = False
    
    def set_objective_lofi(self,fLow,x0):
        self.func = fLow
        self.x0 = x0
        self.vcmObjective = False
    
    def set_objective_vcm(self,fLow,fHigh,x0):
        self.func = ScaledFunction(fLow,fHigh)
        self.x0 = x0
        self.vcmObjective = True

    def add_constraint_lofi(self,gLow):
        self.vcmConstr = False
        self.constr.append(gLow)
    
    def add_constraint_vcm(self,gLow,gHigh):
        self.vcmConstr = True
        self.nVcmConstr += 1
        self.constrVcm.append(ScaledFunction(gLow,gHigh))
    
    def _update_delta(self,rho,err):
        for r in rho:
            if r<=self.eta1 or r>=self.eta3:
                d = self.delta * self.c1
            elif self.eta2<r<self.eta3:
                d = self.delta
            else:
                if err==self.delta:
                    d = self.delta * self.c2
                else:
                    d = self.delta
            if d<self.delta:
                self.delta = d

    def solve(self):
        err = self.tolX + 1.0
        x = self.x0
        nIter = 0
        fnew = 0
        while err>self.tolX or nIter<self.iterMax:
            #print nIter, x, self.delta, fnew
            if self.vcmObjective:
                self.func.construct_scaling_model(x)
            if self.vcmConstr:
                for i in range(self.nVcmConstr):
                    self.constrVcm[i].construct_scaling_model(x)
            tmpCnstrList = self.constr + self.constrVcm
            cnstr = tuple()
            bnds  = tuple()
            for tmpCnstr in tmpCnstrList:
                cnstr = cnstr + ({'type':'ineq','fun':tmpCnstr},)
            if hasattr(x,'__iter__'):
                for xx in x:
                    bnds = bnds + ((xx-self.delta,xx+self.delta),)
            else:
                bnds = ((x-self.delta,x+self.delta),)
            print bnds
            rslt = minimize(self.func,x,method='SLSQP',bounds=bnds,
                            constraints=cnstr)
            xnew = rslt.x
            fnew = rslt.fun
            err = norm(x-xnew)
            x = xnew
            rho = list()
            if self.vcmObjective:
                rho.append(self.func.get_thrust_region_ratio(xnew))
            if self.vcmConstr:
                for i in range(self.nVcmConstr):
                    rho.append(self.constrVcm[i].get_thrust_region_ratio(xnew))
            print rho
            raw_input()
            self._update_delta(rho,err)
            nIter += 1
            