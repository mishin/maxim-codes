# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 21:09:15 2013

@author: Maxim
"""

import sys
#sys.path.append('../')
from vcm import *
from miscTools import normalize, denormalize
from misc_tools import SaveTextData
import airfoil as Af
#from CFDsolver import *
from CFDpolar import *
from FlightConditions import FlightConditions

class AirfoilAnalysis:
    def __init__(self):
        self._clmax = 1.465
        self.clCruise = [0.4,0.6]
        self.thicknessMin = 0.12
        self.thicknessMax = 0.16
        self.af = None
        self.Mcrs = 0.1644
        self.hCrs = 1500.0
        self.Mldg = 0.0726
        self.hLdg = 0.0
        self.cruise = FlightConditions(self.Mcrs, self.hCrs)
        self.landing = FlightConditions(self.Mldg, self.hLdg)
        self.lb = array([0.1, 0.1, 0.1, 0.1, -0.3, -0.3, -0.3])
        self.ub = array([0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1])
        self.zTE = 0.0

    def fLow(self,x):
        self._upd_cst(x)
        self.af.set_trailing_edge(self.zTE)
        try:
            pol = self.af.get_X_polar(self.Mldg,self.landing.Re,alphaSeq=[-10,10,0.25])
            cd = array([pol.get_cd_at_cl(cl) for cl in self.clCruise])
            return cd.mean()
        except:
            self.af.plot()
            return 0.01

    def fLowDeriv(self,x,dx=0.01):
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
        try:
            pol = self.af.get_X_polar(self.Mldg,self.landing.Re,alphaSeq=[0,20,1.0])
            pol.calc_clmax()
            sys.stdout.write('.')
            return -(self._clmax - pol.clmax)
        except:
            sys.stdout.write('x')
            return -0.5

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

    def _upd_cst(self,x,_norm=True):
        if _norm:
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
    def _run_fluent(self,x,_norm=True):
        self._upd_cst(x,_norm)
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
        result = solver.run_for_multiple_aoa(arange(0.,20.,2.0),'ke-realizable')
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


def get_trust_region_ratio(f0, fnew,fHighNew, g0, gnew, gHighNew):
    def get_penalty(f,g):
        P = f
        mu = 1000.
        if hasattr(g,'__iter__'):
            for val in g:
                P += -mu*min(0,val)
        else:
            P = P-mu*min(0,g)
        return P
    rho = (get_penalty(f0,g0) - get_penalty(fHighNew, gHighNew))/(get_penalty(f0,g0) - get_penalty(fnew,gnew))
    return rho

def run_doe_cfd():
    x,f = read_xls_doe()

    aa = AirfoilAnalysis()
    clmax = zeros(len(x))
    clmaxPath = 'KLA100//clmax_DOE3.txt'
    fid = open(clmaxPath,'wt')
    fid.write('Iter\tclmax\talphaClmax\tthickness\n')
    fid.close()
    x0 = array([0.18723832, 0.2479892, 0.26252777, 0.31606257, 0.0819584, -0.11217863, -0.14363534, -0.06480575, -0.27817776, -0.09874038])
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
    c1 = 0.3
    c2 = 2.0
    tol = 1e-6
    gtol = 1e-6
    err = tol+1
    iterMax = 20
    nIter = 0
    delta = 0.4
    # baseline
    #x0 = array([0.18723832, 0.2479892, 0.26252777, 0.31606257, 0.0819584, -0.11217863, -0.14363534, -0.06480575, -0.27817776, -0.09874038])
    x0 = array([0.1670, 0.3340, 0.2601, 0.1971, -0.0713, -0.1951, -0.0498])
    dxu = zeros(len(x0))+0.05
    dxl = zeros(len(x0))-0.05
    #dxu[9] = 0.005
    #dxl[4] = -0.005
    aa.ub = x0+dxu
    aa.lb = x0+dxl
    
    x0 = normalize(x0,aa.lb,aa.ub)
    aa._upd_cst(x0)
    #aa.af.plot()
    gscaled = ScaledFunction(aa.gLow, aa.gHigh,5,scalingType='add')
    #gscaled._initialize_by_doe_points(xdoe,f-aa._clmax)
    
    histFile = 'KLA100//airfoil_root_history_20141206_AVCM.txt'
    fid = open(histFile,'wt')
    fid.write('Iter\tfval\tgHi\tgLow\tgScaled\trho\tdelta\terror\tx\n')
    f0 = aa.fLow(x0)
    g0 = aa.gLow(x0)
    g0h = aa.gHigh(x0)
    out = '%d\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2e\t%.2e\t'%(0,f0,g0+aa._clmax,g0h+aa._clmax,g0h+aa._clmax,0,delta,err)
    for xx in denormalize(x0,aa.lb,aa.ub):
        out += '%.4f\t'%xx
    fid.write(out)
    fid.write('\n')
    fid.close()
    ifevLow = 0
    while err>tol and nIter<iterMax:
        print '--> iteration started'
        nIter += 1
        gscaled.construct_scaling_model(x0)
        bnds = tuple()
        for _x,_xu,_xl in zip(x0,aa.ub,aa.lb):
            bndU = min([_x+delta, _xu])
            bndL = max([_x-delta, _xl])
            bnds += ((bndL, bndU),)
        cnstr = ({'type':'ineq','fun':gscaled},{'type':'ineq','fun':aa.g2},
                 {'type':'ineq','fun':aa.g3})
        rslt = minimize(aa.fLow,x0,method='SLSQP',bounds=bnds,constraints=cnstr,
                        tol=1e-6,jac=aa.fLowDeriv, options={'disp':False})
        xnew = rslt.x
        fnew = rslt.fun
        ifevLow += rslt.nfev
        g0 = gscaled(x0)
        goptScaled = gscaled(xnew)
        goptHigh = gscaled.funcHi(xnew)
        f0 = aa.fLow(x0)
        rho = get_trust_region_ratio(f0,fnew,fnew,g0,goptScaled,goptHigh)
        if goptHigh>=0.0:
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
        #print nIter, fnew, gHinew, aa._clmax, gscaled.funcLo(xnew), gscaled(xnew), rho, delta, err
        out = '%d\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2e\t%.2e\t'%(nIter,fnew,goptHigh+aa._clmax,gscaled.funcLo(xnew)+aa._clmax,gscaled(xnew)+aa._clmax,rho,delta,err)
        for xx in denormalize(xnew,aa.lb,aa.ub):
            out += '%.4f\t'%xx
        fid.write(out)
        fid.write('\n')
        fid.close()
        print out
        print 'lowfi calls: %d, total: %d'%(rslt.nfev, ifevLow)
        if rho>=0:
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
    x0 = array([0.1670, 0.3340, 0.2601, 0.1971, -0.0713, -0.1951, -0.0498])
    x1 = [0.1839,0.3520,0.2512,0.1835,-0.0493,-0.2190,-0.0225]
    x2 = [0.1729,0.3482,0.2707,0.2045,-0.0774,-0.1563,-0.0498]
    
    rsltFile = SaveTextData('af_sub_20141207.txt')
    
    def print_polar(pol):
        print '-->'
        for i in range(len(pol.alpha)):
            print '%.2f\t%.6f\t%.6f\t%.6f'%(pol.alpha[i],pol.cl[i],pol.cd[i],pol.cm[i])

    aa = AirfoilAnalysis()
    aa._upd_cst(x0,False)
    crd0 = copy(aa.af.coord)
    pol0 = aa._run_fluent(x0,False)
    pol0x = aa.af.get_X_polar(aa.Mldg,aa.landing.Re,alphaSeq=[-10,10,1])
    print_polar(pol0)
    print_polar(pol0x)
    
    aa._upd_cst(x1,False)
    crd1 = copy(aa.af.coord)
    pol1 = aa._run_fluent(x1,False)
    pol1x = aa.af.get_X_polar(aa.Mldg,aa.landing.Re,alphaSeq=[-10,10,1])
    print_polar(pol1)
    print_polar(pol1x)
    
    aa._upd_cst(x2,False)
    crd2 = copy(aa.af.coord)
    pol2 = aa._run_fluent(x2,False)
    pol2x = aa.af.get_X_polar(aa.Mldg,aa.landing.Re,alphaSeq=[-10,10,1])
    
    print_polar(pol0)
    print_polar(pol0x)
    print_polar(pol1)
    print_polar(pol1x)
    print_polar(pol2)    
    print_polar(pol2x)
    
    rsltFile.write_array(crd0[:,0],'x0')
    rsltFile.write_array(crd0[:,1],'y0')
    rsltFile.write_array(crd1[:,0],'x1')
    rsltFile.write_array(crd1[:,1],'y1')
    rsltFile.write_array(crd2[:,0],'x2')
    rsltFile.write_array(crd2[:,1],'y2')
    
    rsltFile.write_array(pol0x.cd,'cd0x')
    rsltFile.write_array(pol1x.cd,'cd1x')
    rsltFile.write_array(pol2x.cd,'cd2x')
    rsltFile.write_array(pol0x.cl,'cl0x')
    rsltFile.write_array(pol1x.cl,'cl1x')
    rsltFile.write_array(pol2x.cl,'cl2x')
    
    rsltFile.write_array(pol0.alpha,'a0')
    rsltFile.write_array(pol1.alpha,'a1')
    rsltFile.write_array(pol2.alpha,'a2')
    rsltFile.write_array(pol0.cl,'cl0')
    rsltFile.write_array(pol1.cl,'cl1')
    rsltFile.write_array(pol2.cl,'cl2')
    rsltFile.close()
    
    plt.figure(1)
    plt.hold(True)
    plt.grid(True)
    plt.plot(crd0[:,0],crd0[:,1],'k--',linewidth=2)
    plt.plot(crd1[:,0],crd1[:,1],'r-',linewidth=2)
    plt.plot(crd2[:,0],crd2[:,1],'b:',linewidth=2)
    plt.legend(['baseline','GVFM','VCM'])
    
    plt.figure(2)
    plt.hold(True)
    plt.grid(True)
    plt.axis([0.002, 0.012, -0.5, 1.0])
    plt.plot(pol0x.cd, pol0x.cl, 'ks--',linewidth=2)
    plt.plot(pol1x.cd, pol1x.cl, 'ro-',linewidth=2)
    plt.plot(pol2x.cd, pol2x.cl, 'b^:',linewidth=2)
    plt.xlabel('Drag coefficient')
    plt.ylabel('Lift coefficient')
    plt.legend(['baseline','GVFM','VCM'],'lower right')

    plt.figure(3)
    plt.hold(True)
    plt.grid(True)
    plt.axis([0,20,0,1.8])
    plt.plot(pol0.alpha, pol0.cl, 'ks--',linewidth=2)
    plt.plot(pol1.alpha, pol1.cl, 'ro-',linewidth=2)
    plt.plot(pol2.alpha, pol2.cl, 'b^:',linewidth=2)
    plt.xlabel('Angle of attack, deg')
    plt.ylabel('Lift coefficient')
    plt.legend(['baseline','GVFM','VCM'],'lower right')
    plt.show()
    
if __name__=="__main__":
    #vcm_airfoil_optimization()
    #run_doe_cfd()
    plot_result_af()