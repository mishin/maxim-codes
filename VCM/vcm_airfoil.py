# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 21:09:15 2013

@author: Maxim
"""

import sys
#sys.path.append('../')
from vcm import *
from miscTools import normalize, denormalize
import airfoil as af

class AirfoilAnalysis:
    def __init__(self):
        self._clmax = 1.50
        self.clCruise = [0.4,0.6]
        self.thicknessMin = 0.14
        self.thicknessMax = 0.16
        self.af = None
        self.Mcrs = 0.18
        self.Recrs = 4.4e6
        self.Mldg = 0.09
        self.Reldg = 2.7e6
    
    def fLow(self,x):
        self._upd_cst(x)
        pol = self.af.get_X_polar(self.Mcrs,self.Recrs,alphaSeq=[-5,10,1.0])
        cd = [pol.get_cd_at_cl(cl) for cl in self.clCruise]
        return sum(cd)/len(self.clCruise)
    
    def fLowDeriv(self,x,dx=1e-3):
        grad = zeros(len(x))
        fval = self.fLow(x)
        for i,xx in enumerate(x):
            X = copy(x)
            X[i] = X[i]+dx
            grad[i] = (self.fLow(X)-fval)/dx
        return grad

    def gLow(self,x):
        self._upd_cst(x)
        pol = self.af.get_X_polar(self.Mldg,self.Reldg,alphaSeq=[0,20,1.0])
        pol.calc_clmax()
        sys.stdout.write('.')
        return -(self._clmax - pol.clmax)
    
    def gHigh(self,x):
        self._upd_cst(x)
        pol = self.af.get_J_polar(self.Mldg, self.Reldg)
        pol.calc_clmax()
        sys.stdout.write('-')
        return -(self._clmax - pol.clmax)
    
    def g2(self,x):
        self._upd_cst(x)
        return -(self.thicknessMin - self.af.thickness)

    def g3(self,x):
        self._upd_cst(x)
        return -(self.af.thickness - self.thicknessMax)
    
    def _upd_cst(self,x):
        x = denormalize(x,self.lb,self.ub)
        Au = x[:4]
        Al = array([-x[0],x[4],x[5],x[6]])
        self.af = af.cst(Au,Al)
        self.thickness = self.af.thickness

def read_xls_doe():
    path = 'afInitDoELHC.xls'
    db = af.dbTools.loadDB(path)
    sh = db.selectByName('Sheet1')
    sh = af.dbTools.readDB(sh)
    return sh.readRange(0,0,20)

def vcm_airfoil_optimization():
    aa = AirfoilAnalysis()
    xdoe = read_xls_doe()
    eta1 = 0.25
    eta2 = 0.75
    eta3 = 1.25
    c1 = 0.5
    c2 = 2.0
    tol = 1e-6
    gtol = 1e-3
    err = tol+1
    iterMax = 20
    nIter = 0
    delta = 0.5
    #lb = array([0.05, 0.001, 0.001, 0.001, -0.4, -0.4, -0.4])
    #ub = array([0.4, 0.4, 0.4, 0.4, 0.1, 0.1, 0.1])
    lb = array([0.1, 0.1, 0.1, 0.1, -0.3, -0.3, -0.3])
    ub = array([0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1])
    aa.lb = lb
    aa.ub = ub
    #x0 = (lb+ub)/2.0
    x0 = array([0.1851,0.2578,0.2748,0.1872,-0.1550,-0.2395,0.0593])
    x0 = normalize(x0,lb,ub)
    aa._upd_cst(x0)
    gscaled = ScaledFunction(aa.gLow, aa.gHigh,scalingType='add')
    gscaled._initialize_by_doe_points(xdoe)
    while err>tol and nIter<iterMax:
        print '--> iteration started'
        nIter += 1
        gscaled.construct_scaling_model(x0)
        bnds = tuple()
        for x in x0:
            bnds += ((x-delta, x+delta),)
        cnstr = ({'type':'ineq','fun':gscaled},{'type':'ineq','fun':aa.g2},
                 {'type':'ineq','fun':aa.g3})
        rslt = minimize(aa.fLow,x0,method='SLSQP',bounds=bnds,constraints=cnstr,
                        tol=1e-10,jac=aa.fLowDeriv, options={'disp':False})
        xnew = rslt.x
        fnew = rslt.fun
        gHinew = gscaled.funcHi(xnew)
        if gHinew<=gtol:
            rho = gscaled.get_thrust_region_ratio(xnew,gHinew)
            err = norm(x0-xnew)
        else:
            rho = 1.0
            err = 0.0
        if rho<=eta1 or rho>=eta3:
            delta *= c1
        elif eta1<rho<eta2:
            delta = delta
        else:
            if err==delta:
                delta *= c2
            else:
                delta = delta
        #print nIter,'\t','%.6f\t'%err,fnew,rho,delta
        print '\n%d\t%.2e\t%.6f\t%.2f\t%.4e'%(nIter,err,fnew,rho,delta)
        _xprint = ''
        for xx in xnew: _xprint += '%.4f\t'%xx
        print _xprint
        print gscaled.funcHi(xnew)
        #print gscaled.xPrev,len(gscaled.xPrev)
        x0 = xnew
    
    aa._upd_cst(xnew)
    gscaled.funcHi.display()
    gscaled.funcLo.display()
    print aa.af.thickness
    print aa.gHigh(xnew) + aa._clmax
    pol1 = aa.af.get_J_polar(0.09,2.7e6)
    pol2 = aa.af.get_X_polar(0.17,4.4e6)
    pol1.display()
    pol2.display()
    aa.af.plot()

def plot_result_af():
    x0 = [0.1851,0.2578,0.2748,0.1872,-0.1550,-0.2395,0.0593]
    x1 = [0.1842,0.6701,0.7513,-0.0716,-0.2611,-0.6004,0.8769]
    x2 = [0.2337,0.7288,0.7400,0.0154,-0.2271,-0.4580,1.0059]
    
    aa = AirfoilAnalysis()
    lb = array([0.1, 0.1, 0.1, 0.1, -0.3, -0.3, -0.3])
    ub = array([0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1])
    aa.lb = lb
    aa.ub = ub
    aa._upd_cst(x0)
    crd0 = copy(aa.af.coord)
    aa._upd_cst(x1)
    crd1 = copy(aa.af.coord)
    aa._upd_cst(x2)
    crd2 = copy(aa.af.coord)
    
    plt.figure(1)
    plt.hold(True)
    plt.grid(True)
    #plt.axis('equal')
    plt.plot(crd0[:,0],crd0[:,1],'k--')
    plt.plot(crd1[:,0],crd1[:,1],'r-')
    plt.plot(crd2[:,0],crd2[:,1],'b-')
    plt.legend(['baseline','multiplicative','additive'])
    plt.show()
    
if __name__=="__main__":
    vcm_airfoil_optimization()
    #plot_result_af()