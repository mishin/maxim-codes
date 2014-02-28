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
#from CFDsolver import *
from CFDpolar import *
from FlightConditions import FlightConditions

class AirfoilAnalysis:
    def __init__(self):
        self._clmax = 1.50
        self.clCruise = [0.2,0.3]
        self.thicknessMin = 0.135
        self.thicknessMax = 0.145
        self.af = None
        self.Mcrs = 0.1644
        self.hCrs = 1500.0
        self.Mldg = 0.0726
        self.hLdg = 0.0
        self.cruise = FlightConditions(self.Mcrs, self.hCrs)
        self.landing = FlightConditions(self.Mldg, self.hLdg)
        self.lb = array([0.1, 0.1, 0.1, 0.1, -0.3, -0.3, -0.3])
        self.ub = array([0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1])
        self.zTE = 0.01

    def fLow(self,x):
        self._upd_cst(x)
        self.af.set_trailing_edge(self.zTE)
        try:
            pol = self.af.get_J_polar(self.Mldg,self.landing.Re,alphaSeq=[-10,10,1.0])
            cd = array([pol.get_cd_at_cl(cl) for cl in self.clCruise])
        except ValueError:
            pol.display()
            self.af.plot()
        return cd.mean()

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
        self.af.set_trailing_edge(self.zTE)
        pol = self.af.get_J_polar(self.Mldg,self.landing.Re,alphaSeq=[0,20,1.0])
        pol.calc_clmax()
        sys.stdout.write('.')
        return -(self._clmax - pol.clmax)

    def gHigh(self,x):
        self._upd_cst(x)
        pol = self._run_fluent(x)
        sys.stdout.write('-')
        return -(self._clmax - pol.clmax)

    def g2(self,x):
        self._upd_cst(x)
        return -(self.thicknessMin - self.thickness)

    def g3(self,x):
        self._upd_cst(x)
        return -(self.thickness - self.thicknessMax)

    def _upd_cst(self,x):
        x = denormalize(x,self.lb,self.ub)
        n = len(x)
        if n%2==0:
            Au = x[:n/2]
            Al = x[n/2:]
        else:
            Au = x[:int(n/2)+1]
            Al = hstack([-x[0],x[int(n/2)+1:]])
        self.af = Af.cst(Au,Al,25,'cos')
        self.thickness = self.af.thickness
#    def _run_cfd(self,x):
#        self._upd_cst(x)
#        alphaSeq = array([10.,12,14,16,18])
#        path = paths.CFD_paths()
#        landing = Flight_conditions(0.0,30.0)
#        V = landing.ISA.soundSpeed*0.73
#        landing = Flight_conditions(0.0,V)
#        self.af.create_af_CAT(save=path.file_igs)
#        Airfoil_mesh(path,landing)
#        Airfoil_mesh.yplus_wall = 1.0
#        fluent = Solver(path,landing)
#        fluent.turb_model = 'ke-realizable'
#        for alpha in alphaSeq:
#            fluent.run_fluent(alpha)
#        print fluent.alpha
#        print fluent.cl
#        print fluent.cd
#        print fluent.cm
#        self.af.polar = Af.AirfoilPolar()
#        self.af.polar.Re = landing.Re
#        self.af.polar.Mach = landing.Mach
#        self.af.polar.alpha = fluent.alpha
#        self.af.polar.cl = fluent.cl
#        self.af.polar.cd = fluent.cd
#        self.af.polar.cm = fluent.cm
#        self.af.polar.calc_clmax()
#        path.clean()
#        return self.af.polar
    def _run_fluent(self,x):
        self._upd_cst(x)
        solver = CFDsolver(self.af,self.landing,1.0,mesh='O')
        solver.fluent.residuals['energy']=1e-5
        solver.fluent.relaxationFactor['xvelocity'] = 1e-4
        solver.fluent.residuals['continuity']=1e-4
        solver.mesh._airfoilPts = 75
        solver.mesh._interiorPts = 100
        solver.mesh._dsTE = 2e-4
        solver.mesh._dsLE = 1e-3
        solver.mesh._growthRate = 1.15
        solver.create_mesh()
        result = solver.run_for_multiple_aoa(arange(10.,20.,2.0),'ke-realizable')
        result.Mach = self.landing.Mach
        result.Re = self.landing.Re
        result._calc_clmax()
        return result

def read_xls_doe():
    path = 'KLA_LHC_samples.xls'
    db = Af.dbTools.loadDB(path)
    sh = db.selectByName('Sheet1')
    sh = Af.dbTools.readDB(sh)
    x = sh.readRange(0,0,50)
    return x[:,:-1],x[:,-1]

def run_doe_cfd():
    x,f = read_xls_doe()

    aa = AirfoilAnalysis()
    clmax = zeros(len(x))
    clmaxPath = 'KLA100//clmax_DOE2.txt'
    fid = open(clmaxPath,'wt')
    fid.write('Iter\tclmax\talphaClmax\tthickness\n')
    fid.close()
    x0 = array([0.18723832, 0.2479892, 0.26252777, 0.31606257, 0.0819584, -0.11217863, -0.14363534, -0.06480575, -0.27817776, 0.02874038])
    dxu = zeros(len(x0))+0.05
    dxl = zeros(len(x0))-0.05
    dxu[9] = 0.005
    dxl[4] = -0.005
    aa.ub = x0+dxu
    aa.lb = x0+dxl

    for i,xx in enumerate(x):
        rslt = aa._run_fluent(xx)
        clmax[i] = rslt.clmax
        print '%d\t%.8f\t%.4f\t%.4f\n'%(i,clmax[i],rslt.alphaClmax,aa.thickness)
        fid = open(clmaxPath,'a')
        fid.write('%d\t%.8f\t%.4f\t%.4f\n'%(i+1,clmax[i],rslt.alphaClmax,aa.thickness))
        fid.close()

def vcm_airfoil_optimization():
    aa = AirfoilAnalysis()
    xdoe,f = read_xls_doe()

    eta1 = 0.25
    eta2 = 0.75
    eta3 = 1.25
    c1 = 0.5
    c2 = 2.0
    tol = 1e-6
    gtol = 1e-6
    err = tol+1
    iterMax = 20
    nIter = 0
    delta = 0.4
    # baseline
    x0 = array([0.18723832, 0.2479892, 0.26252777, 0.31606257, 0.0819584, -0.11217863, -0.14363534, -0.06480575, -0.27817776, 0.02874038])
    dxu = zeros(len(x0))+0.05
    dxl = zeros(len(x0))-0.05
    dxu[9] = 0.005
    dxl[4] = -0.005
    aa.ub = x0+dxu
    aa.lb = x0+dxl
    
    x0 = normalize(x0,aa.lb,aa.ub)
    aa._upd_cst(x0)
    gscaled = ScaledFunction(aa.gLow, aa.gHigh,scalingType='add')
    gscaled._initialize_by_doe_points(xdoe,f-aa._clmax)
    
    histFile = 'KLA100//airfoil_root_history_20140227.txt'
    fid = open(histFile,'wt')
    fid.write('Iter\tfval\tgHi\tgLow\tgScaled\trho\tdelta\terror\tx\n')
    fid.close()
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
        rho = gscaled.get_trust_region_ratio(xnew,gHinew)
        if gHinew>=0.0:
            err = 0.0
        else:
            err = norm(x0-xnew)
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
        fid = open(histFile,'a')
#        out = '%d\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2e\t%.2e\t'%(nIter,fnew,gHinew+aa._clmax,gscaled.funcLo(xnew)+aa._clmax,gscaled(xnew)+aa._clmax,rho,delta,err)
#        for xx in denormalize(xnew,aa.lb,aa.ub):
#            out += '%.4f\t'%xx
        #fid.write(out)
        #fid.write('\n')
        fid.close()
        #print out
        x0 = xnew
    
    aa._upd_cst(xnew)
    gscaled.funcHi.display()
    gscaled.funcLo.display()
    print rslt
    print aa.af.thickness
    print aa.gHigh(xnew) + aa._clmax
    #pol1 = aa.af.get_J_polar(0.09,2.7e6)
    pol2 = aa.af.get_X_polar(0.17,4.4e6)
    #pol1.display()
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
    #run_doe_cfd()
    #plot_result_af()