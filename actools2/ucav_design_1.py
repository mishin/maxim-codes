# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 20:58:10 2014

@author: Maxim
"""

import design
from scipy.optimize import fmin_slsqp
import numpy as np
from mission import run_mission_B11, run_mission_B15
from performance import SteadyLevelFlight, ClimbDescent
import copy
from weight_tools import MassComponent

class DesignFormulation(design.Design):
    def setup(self):
        self.lb = np.array([40., 6.0,  0.5,    0.15,   8.0,   2.089, -3.00, 40])
        self.ub = np.array([60., 7.5,  0.7,    0.35,   10.0,  3.600,  1.00, 60])
        self.x0 = np.array([55., 6.91, 0.6006, 0.2656, 9.110, 2.886,  0.00, 55])
        self.xCurrent = np.zeros(len(self.x0))
        self.bnds = np.array([[l,u] for l,u in zip(self.lb,self.ub)])
        # --- constraints ---
        self.WemptyMax      = 3500.0
        self.CnbMin         = 0.0
        self.ClbMax         = -0.05
        self.SMmin          = -0.10
        self.SMmax          = 0.10
        self.RangeMin       = 3000.0
        self.RangeCombatMin = 800.0
        self.RCmin          = 110.0
        self.VmaxMin        = 0.85 # Mach
        self.VminMax        = 0.2
        # --- payload ---
        self.SDB = MassComponent('drop payload', 1132.0, np.array([4.5, 0.0, 0.12]))

    def set_x(self,x):
        """
        To minimize function evaluation, this function re-calculates analysis if 
        new x is given, otherwise precalculated results are used.
        """
        if not (x==self.xCurrent).all():
            self.xCurrent = x
            self._upd_configuration(x)
            try:
                self._upd_analysis()
            except ValueError:
                print x
                print self.g(self.xCurrent)
                print self._cnstrData
                raise ValueError
                
    
    def _upd_configuration(self,x):
        sweepLE = x[0] # leading edge sweep
        cr      = x[1] # root chord
        TR1     = x[2] # 1st segment taper ratio
        TR2     = x[3] # 2nd segment taper ratio
        l       = x[4] # wing span
        l1      = x[5] # 1st segment span (central part)
        twist   = x[6] # tip section twist
        sweep2  = x[7]
        l1 = l1/2.0
        l2 = l/2.0 - l1
        self.set_spans([l1,l2])
        self.set_sweep_angles([sweepLE, sweep2])
        self.set_chord_by_index(cr,0)
        self.set_taper_ratio_by_index(TR1,0)
        self.set_taper_ratio_by_index(TR2,1)
        self.set_twist_by_index(twist,1)
    
    def _upd_analysis(self):
        # update mass and balance
        # update drag
        self._update_mass()
        self._upd_drag()
        V = self.designGoals.cruiseSpeed
        alt = self.designGoals.cruiseAltitude
        self.aero = self.get_aero_trim(V,alt)
        # mission
        ac1 = copy.deepcopy(self)
        ac2 = copy.deepcopy(self)
        if not ac1.mass.payload.item_exists(self.SDB.name):
            ac1.mass.payload.add_component(self.SDB)
        if ac2.mass.payload.item_exists(self.SDB.name):
            ac2.mass.payload.remove_item(self.SDB.name)
        ac3 = copy.copy(ac1)
        ac4 = copy.copy(ac2)
        #print ac1.mass(), ac2.mass(), ac3.mass(), ac4.mass()
        self.R = run_mission_B11(ac2)
        self.Rcombat = run_mission_B15(ac1)
        # performance
        slf = SteadyLevelFlight(ac3)
        clm = ClimbDescent(ac4)
        self.Vmax = slf.run_max_TAS(self.designGoals.cruiseAltitude).Mach
        self.Vmin = slf.run_min_TAS(self.designGoals.cruiseAltitude).Mach
        self.RC = clm.run_max_climb_rate(0).climbRate
        self._cnstrData = np.zeros(10)
        self._cnstrData[0] = self.mass.empty()
        self._cnstrData[1] = self.aero.derivs.Cnb
        self._cnstrData[2] = self.aero.derivs.Clb
        self._cnstrData[3] = self.aero.SM
        self._cnstrData[4] = self.aero.SM
        self._cnstrData[5] = self.R
        self._cnstrData[6] = self.Rcombat
        self._cnstrData[7] = self.RC
        self._cnstrData[8] = self.Vmax
        self._cnstrData[9] = self.Vmin

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
        g = np.zeros(9)
        g[0] = self.WemptyMax - self.mass.empty()
        g[1] = self.aero.derivs.Cnb - self.CnbMin
        g[2] = self.ClbMax - self.aero.derivs.Clb
        g[3] = self.SMmax - self.aero.SM
        g[4] = self.aero.SM - self.SMmin
        g[5] = self.R - self.RangeMin
        g[6] = self.Rcombat - self.RangeCombatMin
        g[7] = self.RC - self.RCmin
        g[8] = self.Vmax - self.VmaxMin
        #g[9] = self.VminMax - self.Vmin
        g[1] *= 1.0e4
        g[2] *= 1.0e2
        g[3] *= 1.0e1
        g[4] *= 1.0e1
        return g
    
    def _set_dvar(self,x):
        sweep1 = x[0]
        sweep2 = x[1]
        cr     = x[2]
        TR1    = x[3]
        TR2    = x[4]
        l1     = x[5]
        l2     = x[6]
        twist1 = x[7]
        twist2 = x[8]
        c2 = cr*TR1
        c3 = c2*TR2
        self.set_chord_by_index(cr,0)
        self.set_chord_by_index(c2,1)
        self.set_chord_by_index(c3,2)
        self.set_spans([l1,l2])
        self.set_sweep_angles([sweep1,sweep2])
        self.set_twist_by_index(twist1,0)
        self.set_twist_by_index(twist2,1)
    
    def run_full_analysis(self,x):
        self._set_dvar(x)
        self._upd_analysis()
        g = np.zeros(8)
        g[0] = self.mass.empty()
        g[1] = self.aero.derivs.Cnb
        g[2] = self.aero.derivs.Clb
        g[3] = self.aero.SM
        g[4] = self.R
        g[5] = self.Rcombat
        g[6] = self.RC
        g[7] = self.Vmax
        return g


def run_optimization():
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    
    rslt = fmin_slsqp(ac.f, ac.x0, f_ieqcons=ac.g, bounds=ac.bnds,
                      epsilon=1e-4,iprint=2)
    
    print ac.g(rslt)
    print ac.aero.display()
    print rslt
    ac.set_x(rslt)
    print ac._cnstrData
    ac.display()

def function_for_sensitivity():
    lb = np.array([40, 40, 6, 0.5, 0.15, 1, 3, -4, -4])
    ub = np.array([60, 60, 7.5, 0.7, 0.35, 1.8, 3.2, 0, 0])
    x0 = np.array([55, 55, 6.91, 0.6006, 0.2656, 1.443, 3.112, 0, -3])
    
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    
    ac._set_dvar(lb)
    ac.display()
    ac._set_dvar(ub)
    ac.display()
    ac._set_dvar(x0)
    ac.display()


if __name__=="__main__":
    function_for_sensitivity()