# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 20:58:10 2014

@author: Maxim
"""

import design
from scipy.optimize import minimize, fmin_slsqp
import numpy as np

class DesignFormulation(design.Design):
    def setup(self):
        self.lb = np.array([40., 6.0,  0.5,    0.15,   8.0,   2.089])
        self.ub = np.array([60., 7.5,  0.7,    0.35,   10.0,  3.600])
        self.x0 = np.array([55., 6.91, 0.6006, 0.2656, 9.110, 2.886])
        self.xCurrent = np.zeros(len(self.x0))
        self.bnds = np.array([[l,u] for l,u in zip(self.lb,self.ub)])
        # --- finite difference ---
        self.dx = 1.0e-4
        # --- constraints ---
        self.WemptyMax = 3533.0
        self.CnbMin = 0.0001
        self.ClbMax = -0.04
        self.SMmin = -0.10
        self.SMmax = 0.10
        self.RangeMin = 3000.0
        self.RangeCombatMin = 800.0
        self.RCmin = 100.0
        self.VmaxMin = 0.7 # Mach
        self.VminMax = 0.2

    def set_x(self,x):
        """
        To minimize function evaluation, this function re-calculates analysis if 
        new x is given, otherwise precalculated results are used.
        """
        if not (x==self.xCurrent).all():
            self.xCurrent = x
            self._upd_configuration(x)
            self._upd_analysis(x)
    
    def _upd_configuration(self,x):
        sweepLE = x[0] # leading edge sweep
        cr      = x[1] # root chord
        TR1     = x[2] # 1st segment taper ratio
        TR2     = x[3] # 2nd segment taper ratio
        l       = x[4] # wing span
        l1      = x[5] # 1st segment span (central part)
        l1 = l1/2.0
        l2 = l/2.0 - l1
        self.set_spans([l1,l2])
        self.set_sweep_angles([sweepLE, sweepLE])
        self.set_chord_by_index(cr,0)
        self.set_taper_ratio_by_index(TR1,0)
        self.set_taper_ratio_by_index(TR2,1)
    
    def _upd_analysis(self,x):
        # update mass and balance
        # update drag
        self._update_mass()
        self._upd_drag()
        V = self.designGoals.cruiseSpeed
        alt = self.designGoals.cruiseAltitude
        self.aero = self.get_aero_trim(V,alt)
    
    def f(self,x):
        self.set_x(x)
        LD = self.aero.coef.CL/self.aero.coef.CD
        return -LD
    
    def _deriv(self,x,dx,func):
        deriv = np.zeros(len(x))
        fval = func(x)
        for i in range(len(x)):
            _x = np.copy(x)
            _x[i] = _x[i]+dx
            deriv[i] = (func(_x)-fval)/dx
        return deriv
    
    def g(self,x):
        self.set_x(x)
        g = np.zeros(5)
        g[0] = self.WemptyMax - self.mass.empty()
        g[1] = self.aero.derivs.Cnb - self.CnbMin
        g[2] = self.ClbMax - self.aero.derivs.Clb
        g[3] = self.SMmax - self.aero.SM
        g[4] = self.aero.SM - self.SMmin
        g = g*1.0e4
        return g


def run_optimization():
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    
    rslt = fmin_slsqp(ac.f, ac.x0, f_ieqcons=ac.g, bounds=ac.bnds,
                      epsilon=5e-4,iprint=2)
    
    print ac.g(rslt)
    print ac.mass.empty()
    print ac.aero.display()
    ac.set_x(rslt)
    ac.display()


if __name__=="__main__":
    run_optimization()