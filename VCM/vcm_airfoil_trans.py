# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 21:09:15 2013

@author: Maxim
"""

import sys
#sys.path.append('../')
from vcm import *
from miscTools import normalize, denormalize
import airfoil as Af
from CFDsolver import *
from tmp_doe import read_samples
import numpy as np

class AirfoilAnalysis:
    def __init__(self):
        self._clmax = 1.50
        self.alphaCruise = 2.0
        self.thicknessMin = 0.07
        self.thicknessMax = 0.13
        self.af = None
        self.Mcrs = 0.20
        self.Recrs = 3e6
        self.lb = array([0.1, 0.1, 0.1, 0.1, -0.3, -0.3, -0.3])
        self.ub = array([0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1])
        self.CmMin = -0.03
    
    def fLow(self,x):
        self._upd_cst(x)
        pol = self.af.get_X_polar(self.Mcrs,self.Recrs,[1,5,1],nIter=100)
        sys.stdout.write('.')
        return -pol.cl[1]/pol.cd[1]
    
    def fHigh(self,x):
        self._upd_cst(x)
        pol = self.af.get_J_polar(self.Mcrs,self.Recrs,[1,5,1])
        sys.stdout.write('*')
        return -pol.cl[1]/pol.cd[1]
    
    def fLowDeriv(self,x,dx=1e-4):
        grad = zeros(len(x))
        fval = self.fLow(x)
        for i,xx in enumerate(x):
            X = copy(x)
            X[i] = X[i]+dx
            grad[i] = (self.fLow(X)-fval)/dx
        return grad

    def gLow(self,x):
        self._upd_cst(x)
        pol = self.af.get_X_polar(self.Mcrs,self.Recrs,[1,5,1],nIter=100)
        sys.stdout.write('-')
        return pol.cm[1]-self.CmMin
    
    def gHigh(self,x):
        self._upd_cst(x)
        pol = self.af.get_J_polar(self.Mcrs,self.Recrs,[1,5,1])
        sys.stdout.write('^')
        return pol.cm[1]-self.CmMin
    def g2(self,x):
        self._upd_cst(x)
        return self.thicknessMax - self.af.thickness
    
    def g3(self,x):
        self._upd_cst(x)
        return self.af.thickness - self.thicknessMin

    def _upd_cst(self,x):
        #x = denormalize(x,self.lb,self.ub)
        Au = x[:4]
        Al = array([-x[0],x[4],x[5],x[6]])
        self.af = Af.cst(Au,Al)
        self.thickness = self.af.thickness

def get_bounds(x0,delta,lb,ub):
    bnds = list()
    for i,xx in enumerate(x0):
        bnd = [max([lb[i],xx-delta]),min([ub[i],xx+delta])]
        bnds.append(bnd)
    return array(bnds,dtype=float)
    
def read_xls_doe():
    path = 'afInitDoELHC.xls'
    db = Af.dbTools.loadDB(path)
    sh = db.selectByName('Sheet1')
    sh = Af.dbTools.readDB(sh)
    x = sh.readRange(0,0,21)
    return x[:,:7],x[:,7]

def transonic_airfoil():
    #x0 = array([0.119087477, 0.160950359,0.203634413,0.192468212,-0.200580639, -0.126010045, 0.107256400e-18])
    x0 = array([0.17042532,0.14831629,0.14576823,0.134351,-0.15162484,-0.13875406,-0.14055989])
    tol  = 1.0e-4
    gtol = 1.0e-4
    maxIter = 50
    xDoe = read_samples('LHC_transonic_af.txt')
    aa = AirfoilAnalysis()
    aa._upd_cst(x0)
    lb = x0 - 0.075
    ub = x0 + 0.075
    err = tol+1.
    niter = 0
    gConverged = False
    xConverged = False
    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.3, 2.0)

    xDoe = denormalize(xDoe,lb,ub,0)
    fscaled = ScaledFunction(aa.fLow, aa.fHigh,0,'add')
    gscaled = ScaledFunction(aa.gLow, aa.gHigh,0,'add')

    fscaled._initialize_by_doe_points(xDoe)
    #gscaled._initialize_by_doe_points(xDoe)
    def print_header():
        print '\nx1\tx2\tf\trho\tdelta\terr\tgScaled\tgHigh'
    while err>tol:
        fscaled.construct_scaling_model(x0)
        #gscaled.construct_scaling_model(x0)
        bnds = get_bounds(x0,delta,lb,ub)
        cnstr = ({'type':'ineq','fun':aa.g3},{'type':'ineq','fun':aa.g2},)
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,constraints=cnstr,tol=1e-10,jac=aa.fLowDeriv)
        #rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-10)
        xnew = rslt.x
        fnew = rslt.fun
        #gScaledNew = gscaled(xnew)
        rho1, fHighNew = fscaled.get_trust_region_ratio(xnew)
        #rho2, gHighNew = gscaled.get_trust_region_ratio(xnew)
        #rho = min([rho1,rho2])
        rho = rho1
        err = np.linalg.norm([x0-xnew])

        delta = trustRegion.adjust(rho,err)
        x0 = xnew
        #print '\n%.4f\t%.4f\t%.4f\t'%(x0[0],x0[1],fnew),'%.2f\t%.2e\t'%(rho,delta),'%.2e'%err
        print '\n',fnew, rho, delta, err
        niter += 1

#        if (-gtol<=gHighNew<=gtol and gScaledNew<=gtol) or (gHighNew>=0.0 and gScaledNew>gtol):
#            gConverged = True
#        else:
#            gConverged = False
    fscaled.funcHi.display()
    fscaled.funcLo.display()
    gscaled.funcHi.display()
    gscaled.funcLo.display()


if __name__=="__main__":
    transonic_airfoil()
