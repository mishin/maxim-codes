# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 13:58:24 2013

@author: Maxim
"""

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