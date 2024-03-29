# -*- coding: utf-8 -*-
"""
Created on Sat Nov 09 13:16:47 2013

@author: Maxim
"""
#import numpy as np
import airfoil as Af
from tmp_doe import read_samples
from vcm import *
import numpy as np
from CFDpolar import *

class AirfoilAnalysis:
    def __init__(self):
        self.alphaCruise = 2.00     
        self.MachCruise = 0.73
        self.thicknessMin = 0.10
        self.thicknessMax = 0.14
        self.af = None
        fc = FlightConditions(0.0,9e3)
        V = fc.atmosphere.soundSpeed * self.MachCruise
        self.fc = FlightConditions(V,9e3)
        self.CmMax = 0.15
        self.hifiHistoryFilePath = 'transonic_af_cfd_history_high_1.txt'
        self.lofiHistoryFilePath = 'transonic_af_cfd_history_low_1.txt'
    
    def _read_cfd_history_file(self,path):
        fid = open(path,'rt')
        lines = fid.readlines()
        fid.close()
        n = len(lines)
        x = zeros([n,7])
        cl = zeros(n)
        cd = zeros(n)
        cm = zeros(n)
        ld = zeros(n)
        for i,line in enumerate(lines):
            if line.strip()!='':
                sLine = line.split()
                x[i,0] = float(sLine[0])
                x[i,1] = float(sLine[1])
                x[i,2] = float(sLine[2])
                x[i,3] = float(sLine[3])
                x[i,4] = float(sLine[4])
                x[i,5] = float(sLine[5])
                x[i,6] = float(sLine[6])
                cl[i] = float(sLine[7])
                cd[i] = float(sLine[8])
                cm[i] = float(sLine[9])
                ld[i] = float(sLine[10])
        return x,cl,cd,cm,ld
    
    def _append_cfd_history_file(self,path,x,cl,cd,cm,ld):
        fid = open(path,'a')
        for _x in x:
            fid.write('%.8f\t'%_x)
        fid.write('%.8f\t%.8f\t%.8f\t%.8f\n'%(cl,cd,cm,ld))
        fid.close()
    
    def _check_cfd_history(self,path,xnew,xtol=1e-8):
        x,cl,cd,cm,ld = self._read_cfd_history_file(path)
        for i,_x in enumerate(x):
            err = np.linalg.norm([xnew-_x])
            if err<=xtol:
                idx = i
                break
        else:
            idx = -1
        return idx
    
    def fLow(self,x):
        path = self.lofiHistoryFilePath
        idx = self._check_cfd_history(path,x)
        if idx==-1:
            cl,cd,cm,ld = self._run_low_fidelity_cfd(x)
            self._append_cfd_history_file(path,x,cl,cd,cm,ld)
            return -ld
        else:
            x,cl,cd,cm,ld = self._read_cfd_history_file(path)
            return -ld[idx]
    
    def fHigh(self,x):
        path = self.hifiHistoryFilePath
        idx = self._check_cfd_history(path,x)
        if idx==-1:
            cl,cd,cm,ld = self._run_high_fidelity_cfd(x)
            self._append_cfd_history_file(path,x,cl,cd,cm,ld)
            return -ld
        else:
            x,cl,cd,cm,ld = self._read_cfd_history_file(path)
            return -ld[idx]
    
    def g1Low(self,x):
        path = self.lofiHistoryFilePath
        idx = self._check_cfd_history(path,x)
        if idx==-1:
            cl,cd,cm,ld = self._run_low_fidelity_cfd(x)
            self._append_cfd_history_file(path,x,cl,cd,cm,ld)
            return self.CmMax - cm
        else:
            x,cl,cd,cm,ld = self._read_cfd_history_file(path)
            return self.CmMax - cm[idx]
    
    def g1High(self,x):
        path = self.hifiHistoryFilePath
        idx = self._check_cfd_history(path,x)
        if idx==-1:
            cl,cd,cm,ld = self._run_high_fidelity_cfd(x)
            self._append_cfd_history_file(path,x,cl,cd,cm,ld)
            return self.CmMax - cm
        else:
            x,cl,cd,cm,ld = self._read_cfd_history_file(path)
            return self.CmMax - cm[idx]
    
    def _run_high_fidelity_cfd(self,x):
        self._upd_cst(x)
        solver = CFDsolver(self.af,self.fc,5,mesh='O')
        solver.fluent.residuals['energy']=1e-6
        solver.fluent.relaxationFactor['xvelocity'] = 1e-3
        solver.mesh._airfoilPts = 220
        solver.mesh._interiorPts = 70
        solver.mesh._dsTE = 2e-4
        solver.mesh._dsLE = 2e-3
        solver.mesh._growthRate = 1.2
        solver.create_mesh()
        result = solver.run_for_single_aoa(self.alphaCruise,iterMax=10000,
                                           turbulenceModel='SA')
        solver.fluent.paths.clean()
        return result.cl, result.cd, result.cm, result.LD
    
    def _run_low_fidelity_cfd(self,x):
        self._upd_cst(x)
        solver = CFDsolver(self.af,self.fc,50,mesh='O')
        solver.fluent.residuals['energy']=1e-4
        solver.fluent.relaxationFactor['xvelocity'] = 1e-3
        solver.mesh._airfoilPts = 55
        solver.mesh._interiorPts = 40
        solver.mesh._dsTE = 5e-4
        solver.mesh._dsLE = 2e-3
        solver.mesh._growthRate = 1.3
        solver.create_mesh()
        result = solver.run_for_single_aoa(self.alphaCruise,iterMax=5000,
                                           turbulenceModel='SA')
        solver.fluent.paths.clean()
        return result.cl, result.cd, result.cm, result.LD
    
    def g2(self,x):
        self._upd_cst(x)
        return self.thicknessMax - self.thickness
    
    def g3(self,x):
        self._upd_cst(x)
        return self.thickness - self.thicknessMin
    
    def _upd_cst(self,x):
        Au = x[:4]
        Al = array([-x[0],x[4],x[5],x[6]])
        self.af = Af.cst(Au,Al)
        self.thickness = self.af.thickness

def read_cfd_output(path,n):
    fid = open(path,'rt')
    lines = fid.readlines()
    fid.close()
    cl = zeros(n)
    cd = zeros(n)
    cm = zeros(n)
    LD = zeros(n)
    for i,line in enumerate(lines[1:]):
        segLine = line.split()
        cl[i] = segLine[0]
        cd[i] = segLine[1]
        cm[i] = segLine[2]
        LD[i] = segLine[3]
    return cl, cd, cm, LD

def get_bounds(x0,delta,lb,ub):
    bnds = list()
    for i,xx in enumerate(x0):
        bnd = [max([lb[i],xx-delta]),min([ub[i],xx+delta])]
        bnds.append(bnd)
    return array(bnds,dtype=float)

    
def get_trust_region_ratio((f0, fnew,fHighNew, g0, gnew, gHighNew)):
    def get_penalty(f,g):
        P = f
        mu = 100.
        if hasattr(g,'__iter__'):
            for val in g:
                P += -mu*min(0,val)
        else:
            P = P-mu*min(0,g)
        return P
    rho = (get_penalty(f0,g0) - get_penalty(fHighNew, gHighNew))/(get_penalty(f0,g0) - get_penalty(fnew,gnew))
    return rho
    
def transonic_airfoil_design():
    x0 = array([1.19087477e-01,1.60950359e-01,2.03634413e-01,1.92468212e-01,-2.00580639e-01,-1.26010045e-01,1.07256400e-18])
    lbDoe = x0 - 0.075
    ubDoe = x0 + 0.075
    lb = x0 - 0.05
    ub = x0 + 0.05
    histPath = 'trans_af_20141107.txt'
    xDoe = read_samples('LHC_transonic_af.txt')
    xDoe = (xDoe+1.0)/2.0*(ubDoe-lbDoe)+lbDoe
    aa = AirfoilAnalysis()
    err = 1.0e-6
    tol = err+1.0
    gtol = 1.0e-6
    maxIter = 20
    nIter = 0
    gConverged = False
    xConverged = False
    
    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.8, 2.0)
    
    clHigh, cdHigh, cmHigh, LDHigh = read_cfd_output('LHC_transonic_CFD_results_omesh_rae2822.txt',20)
    #TODO: precalculate low fidelity solutions
    clLow, cdLow, cmLow, LDLow = read_cfd_output('LHC_transonic_CFD_results_omesh_rae2822_low.txt',20)

    fscaled = ScaledFunction(aa.fLow, aa.fHigh, 0, 'add')
    gscaled = ScaledFunction(aa.g1Low, aa.g1High, 0, 'add')
    
    fscaled._initialize_by_doe_points(xDoe, -LDHigh, -LDLow)
    gscaled._initialize_by_doe_points(xDoe, aa.CmMax-cmHigh, aa.CmMax-cmLow)
    fscaled.dx = 1e-3
    gscaled.dx = 1e-3
    fid = open(histPath,'wt')
    for i in range(len(x0)):
        fid.write('x%d\t'%(i+1))
    fid.write('rno1\trho2\trho\tdelta\terr\tfhighNew\tghighNew\tfnew\tgnew\n')
    fid.close()
    #while xConverged==False or gConverged==False:
    itr = 0
    itrMax = 7
    while itr<itrMax:
        itr += 1
        fscaled.construct_scaling_model(x0)
        gscaled.construct_scaling_model(x0)
        f0 = fscaled(x0)
        g0 = gscaled(x0)
        bnds = get_bounds(x0,delta,lb,ub)
        cnstr = ({'type':'ineq','fun':gscaled,'jac':gscaled.derivative},
                 {'type':'ineq','fun':aa.g2},{'type':'ineq','fun':aa.g3})
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,constraints=cnstr,
                        tol=5e-4,jac=fscaled.derivative)
        xnew = rslt.x
        fnew = rslt.fun
        gnew = gscaled(xnew)
        rho1, fHighNew = fscaled.get_trust_region_ratio(xnew)
        rho2, gHighNew = gscaled.get_trust_region_ratio(xnew)
        #rho = min([rho1,rho2])
        # updated calculation of trust region ratio
        rho = get_trust_region_ratio(f0, fnew,fHighNew, g0, gnew, gHighNew)
        
        err = np.linalg.norm([x0-xnew])
        delta = trustRegion.adjust(rho,err)
        if rho>0.0:
            x0 = xnew
        nIter += 1
        fid = open(histPath,'a')
        for _x in xnew:
            fid.write('%.6f\t'%_x)
        print rho, delta, err, fHighNew, gHighNew
        fid.write('%.4f\t'%rho1)
        fid.write('%.4f\t'%rho2)
        fid.write('%.4f\t'%rho)
        fid.write('%.4e\t'%delta)
        fid.write('%.4e\t'%err)
        fid.write('%.4f\t'%fHighNew)
        fid.write('%.4f\t'%gHighNew)
        fid.write('%.4f\t'%fnew)
        fid.write('%.4f\n'%gnew)
        fid.close()
        gScaledNew = gscaled(xnew)
        if (-gtol<=gHighNew<=gtol and gScaledNew<=gtol) or (gHighNew>=0.0 and gScaledNew>gtol):
            gConverged = True
        else:
            gConverged = False
        
        if nIter>=maxIter or err<=tol:
            xConverged = True
        else:
            xConverged = False
        
    fscaled.funcHi.display()
    fscaled.funcLo.display()
    gscaled.funcHi.display()
    gscaled.funcLo.display()
    
    print xnew
    print fnew
    print fHighNew
    print gHighNew
    print aa.g2(xnew)

if __name__=="__main__":
    transonic_airfoil_design()