# -*- coding: utf-8 -*-
"""
Created on Sat Mar 01 00:12:45 2014

NOTE: optimization is performed on normalized design variables x in [-1;1]

@author: Maxim
"""

from airfoil import *
from misc_tools import Normalization
from flight_conditions import FlightConditions
from CFDpolar import CFDsolver
import sys

class AirfoilObjective:
    def __init__(self,x0,lb,ub):
        self.x0 = x0
        self.af = None
        self.clCruise = [0.2,0.3]
        self.clmaxMin = 1.4
        self.tcMin = 0.135
        self.tcMax = 0.145
        self.norm    = Normalization(lb,ub)
        self.cruise  = FlightConditions(55,1500)
        self.landing = FlightConditions(1.2*20.6, 0)
        self.cdminAlpha = [0,3.0]
        self.afBaseline = cst_x(x0)
    
    def set_cst(self,xnorm):
        x = self.norm.denormalize(xnorm)
#        for x1,x2 in zip(xnorm,x):
#            print '%.4f\t%.4f'%(x1,x2)
#        raw_input()
        self.af = cst_x(x,35)
        self.tc = self.af.thickness
    
    def f(self,x):
        sys.stdout.write('.')
        self.set_cst(x)
        pol = self.af.get_jfoil_polar(self.cruise.Mach, self.cruise.Re, [-10,10,0.5])
        try:
            cd = np.array([pol.get_cd_at_cl(cl) for cl in self.clCruise])
            return cd.mean()
        except ValueError:
            #self.af.display()
            #pol.display()
            return pol.cd[0]
            
        
    
    def fDeriv(self,x,dx=1e-2):
        grad = np.zeros(len(x))
        fval = self.f(x)
        for i,xx in enumerate(x):
            X = np.copy(x)
            X[i] = X[i]+dx
            grad[i] = (self.f(X)-fval)/dx
        return grad
    
    def g1Deriv(self,x,dx=2e-3):
        grad = np.zeros(len(x))
        fval = self.g1high(x)
        for i,xx in enumerate(x):
            X = np.copy(x)
            X[i] = X[i]+dx
            grad[i] = (self.g1high(X)-fval)/dx
        return grad
        
    def run_cfd(self,x,alpha=[12,14,16,18]):
        self.set_cst(x)
        self.af.coord = self.af.pts
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
        result = solver.run_for_multiple_aoa(alpha,'ke-realizable')
        result.Mach = self.landing.Mach
        result.Re = self.landing.Re
        result._calc_clmax()
        sys.stdout.write('|')
        return result
    
    def g1high(self,x):
        pol = self.run_cfd(x)
        return pol.clmax - self.clmaxMin
    
    def g1low(self,x):
        sys.stdout.write('.')
        self.set_cst(x)
        pol = self.af.get_jfoil_polar(self.landing.Mach, self.landing.Re, [0,20,1])
        return pol.clmax-self.clmaxMin

    def g2(self,x):
        self.set_cst(x)
        return self.tc - self.tcMin

    def g3(self,x):
        self.set_cst(x)
        return self.tcMax - self.tc
    
    def g4(self,x):
        self.set_cst(x)
        pol = self.af.get_jfoil_polar(self.cruise.Mach, self.cruise.Re, [-10,10,0.5])
        return pol.get_alpha_cdmin() - self.cdminAlpha[0]
    def g5(self,x):
        self.set_cst(x)
        pol = self.af.get_jfoil_polar(self.cruise.Mach, self.cruise.Re, [-10,10,0.5])
        return self.cdminAlpha[1] - pol.get_alpha_cdmin()


def run_test1():
    x0 = np.array([0.18723832, 0.2479892, 0.26252777, 0.31606257, 0.0819584, -0.11217863, -0.14363534, -0.06480575, -0.27817776, -0.09874038])
    dxu = np.zeros(len(x0))+0.05
    dxl = np.zeros(len(x0))-0.05
    dxu[9] = 0.005
    dxl[4] = -0.005
    ub = x0+dxu
    lb = x0+dxl
    objective = AirfoilObjective(x0,lb,ub)
    x0 = objective.norm.normalize(x0)
    objective.set_cst(x0)
    print objective.cruise.Re
    print objective.landing.Re
    objective.af.display()
    print objective.g1high(x0)
    print objective.f(x0)
    print objective.g1low(x0)
    print objective.g2(x0)
    print objective.g3(x0)


if __name__=="__main__":
    run_test1()