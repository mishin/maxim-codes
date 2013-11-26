# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 12:56:38 2013

@author: Maxim
"""

from functionHandling import *
from scipy.interpolate import Rbf
from scipy.optimize import minimize

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

class ScalingFunction:
    def __init__(self,funcHigh,funcLow,scalingType='add',warmUpRuns=3,onedim=False,dx=1e-4):
        self.type = scalingType
        self.oneDim = onedim
        self._dx = dx
        self.warmUpRuns = warmUpRuns
        self._histScaling = list()
        if onedim:
            self.funcHigh = Function1D(funcHigh)
            self.funcLow = Function1D(funcLow)
        else:
            self.funcHigh = FunctionND(funcHigh)
            self.funcLow = FunctionND(funcLow)
    
    def __call__(self,x,save=True):
        if self.type=='add':
            return self.scalingFactor(x) + self.funcLow(x,save)
        elif self.type=='mult':
            return self.scalingFactor(x) * self.funcLow(x,save)

    def construct_scaling_model(self,x0,fHigh=None,fLow=None):
        self.x0 = x0
        if fHigh==None:
            fHigh = self.funcHigh(x0)
        if fLow==None :
            fLow  = self.funcLow(x0)
        scalingFactor = self._get_scaling_factor(fHigh,fLow)
        self._histScaling = np.array([self._get_scaling_factor(fh,fl) for fh,fl in zip(self.funcHigh._histFpart,self.funcLow._histFpart)])
        if len(self.funcHigh._histFpart)<self.warmUpRuns:
            scalingGradient = self._get_scaling_gradient(x0,fHigh,fLow)
            self.scalingFactor = self._get_taylor_series(x0,scalingFactor,scalingGradient)
        else:
            self.scalingFactor = RbfMod(self.funcHigh._histXpart,self._histScaling)

    def _get_scaling_gradient(self,x0,fHigh,fLow):
        fHigh,gradHigh = self.funcHigh.get_gradient(x0,fval=fHigh)
        fLow, gradLow  = self.funcLow.get_gradient(x0,fval=fLow)
        if self.type=='add':
            return gradHigh-gradLow
        elif self.type=='mult':
            return (gradHigh*fLow - gradLow*fHigh)/(fLow*fLow)
    
    def _get_scaling_factor(self,fHigh,fLow):
        if self.type=='add':
            return fHigh-fLow
        elif self.type=='mult':
            return fHigh/fLow
    
    def _get_taylor_series(self,x0,scalingFactor,scalingGradient):
        if self.oneDim:
            return Taylor1D(x0,scalingFactor,scalingGradient)
        else:
            return TaylorND(x0,scalingFactor,scalingGradient)
    
    def get_trust_region_ratio(self,x):
        fscaled = self.__call__(x)
        fscaled0 = self.__call__(self.x0)
        fHigh = self.funcHigh(x)
        fHigh0 = self.funcHigh(self.x0)
        if fscaled==fscaled0:
            return float('inf')
        else:
            return (fHigh0 - fHigh) / (fscaled0 - fscaled)
    
    def initialize_by_doe(self,xDOE):
        for x in xDOE[:-1]:
            self.funcHigh(x,True)
            self.funcLow(x,True)
        self.construct_scaling_model(xDOE[-1])

    def initialize_from_file(self,filePathHigh,filePathLow=None):
        self.funcHigh.read_history(filePathHigh)
        xDOE = self.funcHigh._histXpart
        if filePathLow==None:
            for x in xDOE:
                self.funcLow(x)
        else:
            self.funcLow.read_history(filePathLow)
        fHigh = self.funcHigh._histFpart[-1]
        fLow = self.funcLow._histFpart[-1]
        self.construct_scaling_model(xDOE[-1],fHigh,fLow)
    
    def get_gradient(self,x):
        if self.oneDim:
            return (self.__call__(x+self._dx) - self.__call__(x)) / self._dx
        else:
            fval = self.__call__(x)
            grad = np.zeros(len(x))
            for i in range(len(x)):
                X = np.copy(x)
                X[i] = X[i]+self._dx
                grad[i] = (self.__call__(X,False)-fval)/self._dx
            return grad


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


class VFM:
    def __init__(self,x0,xL,xU,delta=0.0,tolX=1e-6,tolG=1e-6,iterMax=20,warmUp=3):
        
        self._warmUpRuns = warmUp
        self._histScaled   = list()
        self._histScFactor = list()
        self._histRho      = list()
        self._histDelta    = list()
        self._histFilePath = None
        self._objective    = None
        self._vfmObjective = True
        self._cnstrVfm     = list()
        self._cnstrExact   = list()
        self.x0 = x0
        self.xL = xL
        self.xU = xU
        if hasattr(x0,'__iter__'):
            self._ndim = len(x0)
        else:
            self._ndim = 1
        if self._ndim==1:
            self.onedim = True
        else:
            self.onedim = False
        self._nVfmConstr = 0
        self.isConstrained = False
        self.tolX = tolX
        self.tolG = tolG
        self.tolXopt = self.tolX*1e-2
        self.iterMax = iterMax
        if delta==0.0:
            if self.onedim:
                delta = min([min(self.xU-self.x0,self.x0-self.xL)])
            else:
                delta = min([min(xu-x,x-xl) for x,xu,xl in zip(self.x0,self.xU,self.xL)])
        self.delta0 = delta
        self.trm = TrustRegionManagement(self.delta0)

    def set_objective_exact(self,func):
        self._objective = func
        self._vfmObjective = False

    def set_objective_approximate(self,funcHigh,funcLow,scalingType='add',useJac=False):
        self._objective = ScalingFunction(funcHigh,funcLow,scalingType,self._warmUpRuns,self.onedim)
        self._useJac = useJac
        if useJac:
            self._jac = ScalingFunction.get_gradient
        self._vfmObjective = True

    def add_constraint_exact(self,func):
        self._cnstrExact.append(func)
        self.isConstrained = True

    def add_constraint_approximate(self,funcHigh,funcLow,scalingType='add'):
        g = ScalingFunction(funcHigh,funcLow,scalingType,self._warmUpRuns,self.onedim)
        self._cnstrVfm.append(g)
        self.isConstrained = True
        self._nVfmConstr += 1

    def solve(self):
        if self.isConstrained:
            return self._solve_constrained()
        else:
            return self._solve_unconstrained()
    
    def _solve_unconstrained(self):
        err = self.tolX + 1.0
        nIter = 0
        x0 = self.x0
        delta = self.delta0
        fsc = self._objective
        while err>self.tolX and nIter<self.iterMax:
            fsc.construct_scaling_model(x0)
            bnds = self._get_bounds(x0,delta)
            if self._useJac:
                rslt = minimize(fsc,x0,method='SLSQP',bounds=bnds,
                                tol=self.tolXopt,jac=self._jac)
            else:
                rslt = minimize(fsc,x0,method='SLSQP',bounds=bnds,
                                tol=self.tolXopt)
            xnew = rslt.x
            fnew = rslt.fun
            rho = fsc.get_trust_region_ratio(xnew)
            err = np.linalg.norm([x0-xnew])
            if self._histFilePath!=None:
                self._write_unconstrained_history(nIter,xnew,fnew,fsc.funcHigh(xnew),
                                                  fsc.funcLow(xnew),rho,delta,err)
            delta,x0 = self.trm.adjust(rho,x0,xnew)
            nIter += 1
        return xnew

    def _solve_constrained(self):
        pass
    
    def _get_bounds(self,x,delta):
        bnds = np.zeros([self._ndim,2])
        if not hasattr(x,'__iter__'):
            bnds = np.zeros([1,2])
            bnds[0,0] = max(self.xL, x-delta)
            bnds[0,1] = min(self.xU, x+delta)
        else:
            for i,xx in enumerate(x):
                bnds[i,0] = max(self.xL[i],xx-delta)
                bnds[i,1] = min(self.xU[i],xx+delta)
        return bnds
    
    def _write_unconstrained_history(self,nIter,x,fSc,fHi,fLo,rho,delta,err):
        fid = open(self._histFilePath,'a')
        if nIter==0:
            fid.write('Iter\t')
            for i in range(len(x)):
                fid.write('x%d\t'%i)
            fid.write('fScaled\tfHigh\tfLow\trho\tdelta\terr\n')
        fid.write('%d\t'%nIter)
        for _x in x:
            fid.write('%.6f\t'%_x)
        fid.write('%.6f\t%.6f\t%.6f\t%.4f\t%.4e\t%.4e\n'%(fSc,fHi,fLo,rho,delta,err))
        fid.close()

def run_test1():
    x = np.linspace(0,1,50)
    fhigh = forrester
    flow = lambda x: 0.5*fhigh(x) + 10.*(x-.5)+5.
    
    fsc = ScalingFunction(fhigh,flow,'mult',3,True)
    fsc.construct_scaling_model(0.1)
    fsc.construct_scaling_model(0.9)
    fsc.construct_scaling_model(0.5)
    print fsc.get_trust_region_ratio(0.8)

    y1 = np.array([fhigh(_x) for _x in x])
    y2 = np.array([flow(_x) for _x in x])
    y3 = np.array([fsc(_x) for _x in x])
    plt.figure(1)
    plt.hold(True)
    plt.plot(x,y1,'r-')
    plt.plot(x,y2,'b-')
    plt.plot(x,y3,'g-')
    plt.plot(fsc.funcHigh._histXpart, fsc.funcHigh._histFpart,'ko')
    plt.show()
    

def run_test2():
    fhigh = lambda x: forrester(np.linalg.norm(x))
    flow = lambda x: fhigh(x)+20*np.linalg.norm(x) + 30*sum(x)
    
    fsc = ScalingFunction(fhigh,flow,'add',3,False)
#    fsc.construct_scaling_model(np.array([0.5,0.8]))
#    fsc.construct_scaling_model(np.array([0.8,0.5]))
#    fsc.construct_scaling_model(np.array([0.3,0.4]))
#    fsc.construct_scaling_model(np.array([0.9,0.1]))
    xDOE = np.array([[0.5,0.5],[0.8,0.8],[0.9,0.1]])
    #fsc.initialize_by_doe(xDOE)
    fsc.initialize_from_file('tmp_high_historyZ.txt','tmp_low_historyZ.txt')
    #fsc.funcHigh.write_history('tmp_high_historyZ.txt')
    #fsc.funcLow.write_history('tmp_low_historyZ.txt')
    print fsc.get_trust_region_ratio(np.array([0.45,0.48]))
    lb = [0, 0.0]
    ub = [1., 1.]
    dx = 0.05
    x = y = np.arange(lb[0],ub[0],dx)
    X, Y = np.meshgrid(x,y)
    Z1 = np.array([fsc.funcHigh(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z2 = np.array([fsc.funcLow(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
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
    plt.show()

def run_vfm1():
    x = np.linspace(0,1,50)
    def fhigh(x):
        #x = x[0]
        return (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
    def flow(x):
        #x = x[0]
        return 0.5*fhigh(x) + 10.*(x-.5)+5.
    
    vfm = VFM(0.5,0.0,1.0)
    vfm.set_objective_approximate(fhigh,flow,'add')
    vfm._histFilePath = 'VFM_1dTest_run.txt'
    vfm.solve()

if __name__=="__main__":
    run_vfm1()