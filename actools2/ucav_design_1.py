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
        self.lb = np.array([40, 40, 6.00, 3.00, 0.5, 1.00, 3.000, -4, -4])
        self.ub = np.array([60, 60, 7.50, 5.25, 1.8, 1.80, 3.200,  0,  0])
        self.x0 = np.array([55, 55, 6.91, 4.15, 1.1, 1.44, 3.115,  0, -3])
        self.xCurrent = np.zeros(len(self.x0))
        self.bnds = np.array([[l,u] for l,u in zip(self.lb,self.ub)])
        # --- constraints ---
        self.WemptyMax      = 3400.0
        self.CnbMin         = 0.0001
        self.ClbMax         = -0.05
        self.SMmin          = -0.05
        self.SMmax          = 0.10
        self.combatRadiusMin= 900.0
        self.RCmin          = 125.0
        self.VmaxMin        = 0.90 # Mach
        # --- payload ---
        self.SDB = MassComponent('drop payload', 1132.0, np.array([4.5, 0.0, 0.12]))
        # --- misc ---
        self.analysisData = np.zeros(8)

    def set_x(self,x,fullUpdate=True):
        """
        To minimize function evaluation, this function re-calculates analysis if 
        new x is given, otherwise precalculated results are used.
        """
        if not (x==self.xCurrent).all():
            self.xCurrent = x
            self._upd_configuration(x)
            try:
                self._upd_analysis(fullUpdate)
            except:
                print x
                print self.analysisData
                raise ValueError
    
    def _upd_configuration(self,x):
        sweep1 = x[0]
        sweep2 = x[1]
        chord1 = x[2]
        chord2 = x[3]
        chord3 = x[4]
        span1  = x[5]
        span2  = x[6]
        twist1 = x[7]
        twist2 = x[8]

        self.set_spans([span1,span2])
        self.set_sweep_angles([sweep1, sweep2])
        self.set_chords([chord1,chord2,chord3])
        self.set_twist_by_index(twist1,0)
        self.set_twist_by_index(twist2,1)
    
    def _upd_analysis(self,fullUpdate):
        self._update_mass()
        self._upd_drag()
        V   = self.designGoals.cruiseSpeed
        alt = self.designGoals.cruiseAltitude
        self.aero = self.get_aero_single_point(V,alt,0)
        self.analysisData[0] = self.aero.coef.CL/self.aero.coef.CD

        if fullUpdate:
        # mission
            ac1 = copy.deepcopy(self)
            ac2 = copy.deepcopy(self)
    
            Wf = ac1.mass.fuel.mass
            CGf = ac1.mass.fuel.coords
            if not ac1.mass.payload.item_exists(self.SDB.name):
                ac1.mass.payload.add_component(self.SDB)
            if ac2.mass.payload.item_exists(self.SDB.name):
                ac2.mass.payload.remove_item(self.SDB.name)
    
            self.combatRadius = run_mission_B15(ac1)
            # performance
            ac1.mass.set_fuel_mass(Wf,CGf)
            ac2.mass.set_fuel_mass(Wf,CGf)
            if not ac1.mass.payload.item_exists(self.SDB.name):
                ac1.mass.payload.add_component(self.SDB)
            if not ac2.mass.payload.item_exists(self.SDB.name):
                ac2.mass.payload.add_component(self.SDB)
    
            slf = SteadyLevelFlight(ac1)
            clm = ClimbDescent(ac2)
            
            
            self.analysisData[1] = self.mass.empty()
            self.analysisData[2] = self.aero.derivs.Cnb
            self.analysisData[3] = self.aero.derivs.Clb
            self.analysisData[4] = self.aero.SM
            self.analysisData[5] = self.combatRadius
            self.analysisData[6] = clm.run_max_climb_rate(0).climbRate
            self.analysisData[7] = slf.run_max_TAS(alt).Mach

    def f(self,x):
        self.set_x(x,False)
        return -self.analysisData[0]

    def g(self,x):
        self.set_x(x,True)
        g = np.zeros(8)
        g[0] = self.WemptyMax - self.analysisData[1]
        g[1] = self.analysisData[2] - self.CnbMin
        g[2] = self.ClbMax - self.analysisData[3]
        g[3] = self.SMmax - self.analysisData[4]
        g[4] = self.analysisData[4] - self.SMmin
        g[5] = self.analysisData[5] - self.combatRadiusMin
        g[6] = self.analysisData[6] - self.RCmin
        g[7] = self.analysisData[7] - self.VmaxMin
#
        g[0] *= 1e-2
        g[1] *= 1e4
        g[2] *= 1e3
        g[3] *= 1e3
        g[4] *= 1e3
        g[5] *= 1e-2
        g[6] *= 1e-1
        g[7] *= 1e2
        
        return g


def run_optimization():
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.propulsion._build_thrust_table()
    ac.setup()
    ac.set_x(ac.x0)

    rslt = fmin_slsqp(ac.f, ac.x0, f_ieqcons=ac.g, bounds=ac.bnds,
                      epsilon=1e-3,iprint=2,acc=2e-3)
    ac.set_x(rslt)
    print ac.analysisData
    print ac.g(rslt)
    print ac.aero.display()
    print rslt
    ac.display()


def function_for_sensitivity():
    pathIn = 'ucav_input.txt'
    pathOut = 'ucav_output.txt'

    fid = open(pathIn,'rt')
    line = fid.readline().split()
    fid.close()

    n = len(line)
    x = np.zeros(n)

    for i in range(n):
        x[i] = float(line[i])
#    lb = np.array([40, 40, 6, 0.5, 0.15, 1, 3, -4, -4])
#    ub = np.array([60, 60, 7.5, 0.7, 0.35, 1.8, 3.2, 0, 0])
#    x0 = np.array([55, 55, 6.91, 0.6006, 0.2656, 1.443, 3.112, 0, -3])
    
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    out = ac.run_full_analysis(x)
    print out
    return out


def run_design_table():
    from misc_tools import read_tabulated_data_without_header, Normalization
    data = read_tabulated_data_without_header('DOE_LHS250_FFD2_3.txt')

    pathOut = 'design_out3.txt'
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.propulsion._build_thrust_table()
    ac.setup()

    normalizers = list()
    for l,u in zip(ac.lb, ac.ub):
        normalizers.append(Normalization(l,u))

    dataNew = np.zeros(data.shape)
    n = data.shape[1]
    for i,xNorm in enumerate(data):
        dataNew[i] = np.array([normalizers[ii].denormalize(xNorm[ii]) for ii in range(n)])
        ac.set_x(dataNew[i])
        fid = open(pathOut,'at')
        for val in ac.analysisData:
            fid.write('%.20f\t'%val)
        fid.write('\n')
        fid.close()
        print i,'\n', ac.analysisData


if __name__=="__main__":
#    ac = DesignFormulation()
#    ac.load_xls('Baseline1')
#    ac.propulsion._build_thrust_table()
#    ac.setup()
#    ac.set_x(ac.x0)
    run_design_table()