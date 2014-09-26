# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 20:58:10 2014

@author: Maxim
"""

import design
from scipy.optimize import fmin_slsqp
import numpy as np
from mission import run_mission_B15
from performance import SteadyLevelFlight, ClimbDescent
import copy
from weight_tools import MassComponent

from misc_tools import read_tabulated_data_without_header, Normalization, RbfMod


class DesignFormulation(design.Design):
    def setup(self):
        self.neval = 0
        self.lb = np.array([40, 40, 6.00, 3.00, 0.5, 1.00, 3.000, -4, -4])
        self.ub = np.array([60, 60, 7.50, 5.25, 1.8, 1.80, 3.200,  0,  0])
        self.x0 = np.array([55, 55, 6.91, 4.15, 1.1, 1.44, 3.115,  0, -3])
        #self.x0 = np.array([50.146128, 44.359194, 6.840497, 3.000000, 0.500000, 1.800000, 3.117574, 0.000000, -3.000000])
        self.norm = Normalization(self.lb, self.ub,-1.,1.)
        self.xCurrent = np.zeros(len(self.x0))
        self.x0norm = self.norm.normalize(self.x0)
        #self.bnds = np.array([[l,u] for l,u in zip(self.lb,self.ub)])
        # --- constraints ---
        self.WemptyMax      = 3400.0
        self.CnbMin         = 3e-4
        self.ClbMax         = -6.0e-2
        self.SMmin          = -0.05
        self.SMmax          = 0.10
        self.combatRadiusMin= 750.0
        self.RCmin          = 125.0
        self.VmaxMin        = 0.90 # Mach
        self.minSectionLen  = 1.4*self.propulsion.engine.length
        self.CmdeReq        = -0.008
        # --- payload ---
        self.SDB = MassComponent('drop payload', 1132.0, np.array([4.5, 0.0, 0.12]))
        # --- misc ---
        self.analysisData = np.zeros(10)
        self._upd_approximation()
        # --- initial data ---
        self.Wf  = self.mass.fuel.mass
        self.CGf = self.mass.fuel.coords
        #self.engineOffset = 0.75*self.propulsion.engine.length
        self.engineOffset = 0.75*self.propulsion.engine.length
        self.fileFeasible = 'ucav_feasible.txt'

    def _upd_approximation(self):
        pathIn = 'design_out4.txt'
        pathInSamples = 'DOE_LHS250_FFD2_3_3.txt'
        learnData = read_tabulated_data_without_header(pathIn)
        xNorm = read_tabulated_data_without_header(pathInSamples)
        learnData = np.transpose(learnData)
        _Cnb     = RbfMod(xNorm, (learnData[2])*1e3)
        self.Cnb = lambda x: _Cnb(x)/1e3
        self.Clb = RbfMod(xNorm, learnData[3])
        self.LD  = RbfMod(xNorm, learnData[0])

    def set_x(self,xnorm,fullUpdate=True):
        """
        To minimize function evaluation, this function re-calculates analysis if 
        new x is given, otherwise precalculated results are used.
        """
        x = self.norm.denormalize(xnorm)
        #if not (x==self.xCurrent).all():
        #self.xCurrent = x
        try:
            self._upd_configuration(x)
            self._upd_analysis(xnorm,fullUpdate)
        except:
            print x
            print xnorm
            self._upd_configuration(x)
            print 'Error occured'
            self.display_2d()
            self._upd_analysis(xnorm,fullUpdate)
#                self.analysisData[0] = 0.0
#                self.analysisData[1:] = -np.ones(len(self.analysisData)-1)*100.0
#                raise ValueError
    
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
        self._set_engine_location()
        self._update_mass()
    
    def _set_engine_location(self):
        engineDiameter = self.propulsion.engine.diameter
        x,l = self.wing.get_max_segment_length(engineDiameter)
        self.maxSectionLength = l
        cgNew = x+self.engineOffset
        #self.mass.empty.update_item_cg('engine',cgNew,0,0)
        self.set_engine_cg(cgNew,0,0)
        self.mass.payload.update_item_cg('drop payload',cgNew,.0,.0)
    
    def _upd_analysis(self,xnorm,fullUpdate):
        self._update_mass()
        self._upd_drag()
        V   = self.designGoals.cruiseSpeed
        alt = self.designGoals.cruiseAltitude
        self.aero = self.get_aero_single_point(V,alt,2.0)
        self.analysisData[0] = self.aero.coef.CL/self.aero.coef.CD
        #self.analysisData[0] = self.LD(xnorm) #self.aero.coef.CL/self.aero.coef.CD

        if fullUpdate:
        # mission
            self.mass.set_fuel_mass(self.Wf, self.CGf)
            if not self.mass.payload.item_exists(self.SDB.name):
                self.mass.payload.add_component(self.SDB)
            self.mass.set_fuel_mass(self.Wf,self.CGf)
            self.combatRadius = run_mission_B15(self) /1e3
            # performance
            self.mass.set_fuel_mass(self.Wf,self.CGf)
            if not self.mass.payload.item_exists(self.SDB.name):
                self.mass.payload.add_component(self.SDB)

            slf = SteadyLevelFlight(self)
            self.analysisData[1] = self.mass.empty()
            self.analysisData[2] = self.aero.derivs.Cnb
            self.analysisData[3] = self.aero.derivs.Clb  #self.aero.derivs.Clb
            self.analysisData[4] = self.aero.SM
            self.analysisData[5] = self.combatRadius
            self.analysisData[7] = slf.run_max_TAS(alt).Mach
            S = self.wing.area
            q = self.designGoals.fc.dynamicPressure
            self.analysisData[8] = self.aero.coef.CL*q*S

            self.mass.set_fuel_mass(self.Wf,self.CGf)
            if not self.mass.payload.item_exists(self.SDB.name):
                self.mass.payload.add_component(self.SDB)
            clm = ClimbDescent(self)
            self.analysisData[6] = clm.run_max_climb_rate(0).climbRate
            self.analysisData[9] = self.aero.derivs.Cmde
        print self.analysisData[0]

    def f(self,x):
        self.neval += 1
        self.set_x(x,False)
        return -self.analysisData[0]

    def g(self,x):
        self.set_x(x,True)
        _x = self.norm.denormalize(x)
        g = np.zeros(12)
        g[0] = self.WemptyMax - self.analysisData[1]
        g[1] = self.analysisData[2] - self.CnbMin
        g[2] = self.ClbMax - self.analysisData[3]
        g[3] = self.SMmax - self.analysisData[4]
        g[4] = self.analysisData[4] - self.SMmin
        g[5] = self.analysisData[5] - self.combatRadiusMin
        g[6] = self.analysisData[6] - self.RCmin
        g[7] = self.analysisData[7] - self.VmaxMin
        g[8] = self.maxSectionLength - self.minSectionLen
        g[9] = self.analysisData[8] - self.get_mass()*9.81
        g[10] = _x[0] - _x[1]
        g[11] = self.CmdeReq - self.analysisData[9]
        print g>0
        out = ''
        for v in g:
            out += '%.0e '%v
        print out
        print '%.4f\t%.4f\t%.4f'%(self.analysisData[1],self.analysisData[5],self.analysisData[8])
        
#        if all(g>0):
#            fid = open(self.fileFeasible,'at')
#            for _x in self.norm.denormalize(x):
#                fid.write('%.6f\t'%_x)
#            for val in self.analysisData:
#                fid.write('%.4e\t'%val)
#            fid.write('\n')
#            fid.close()

        g[0] *= 1e-2
        g[1] *= 1e4
        g[2] *= 1e3
        g[3] *= 1e3
        g[4] *= 1e3
        g[5] *= 1e-1
        g[6] *= 1e-1
        g[7] *= 1e2
        g[11] *= 1e3
        return g*100.


def run_optimization():
    import numdifftools as nd
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    #ac.propulsion._build_thrust_table()
    ac.setup()
#    ac.set_x(ac.x0norm)
#    ac.display()
    
    bnds = np.ones([len(ac.x0),2])
    bnds[:,0] = -bnds[:,0]
    ac.set_x(ac.norm.normalize(ac.x0))
    
    fd = nd.Gradient(ac.f,step_max=1e-2, romberg_terms=1)
    gd = nd.Jacobian(ac.g,step_max=1e-2, romberg_terms=1, vectorized=True)

#    print gd(ac.x0norm)
#    print ac.neval
#    print 'completed gradient calculation'
    #raw_input()
#
#    print ac.analysisData
#    #raw_input()
    rslt = fmin_slsqp(ac.f, ac.x0norm, f_ieqcons=ac.g, bounds=bnds, iprint=2,epsilon=1e-2)#,
                      #fprime=fd, fprime_ieqcons=gd)
    ac.set_x(rslt)
    print ac.g(rslt)>=0.0
    print 'LD', ac.analysisData[0]
    print 'We', ac.analysisData[1]
    print 'Cnb',ac.analysisData[2]
    print 'Clb',ac.analysisData[3]
    print 'SM',ac.analysisData[4]
    print 'R',ac.analysisData[5]
    print 'RC',ac.analysisData[6]
    print 'M',ac.analysisData[7]

    ac.aero.display()
    print rslt
    print ac.norm.denormalize(rslt)
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
    data = read_tabulated_data_without_header('DOE_LHS250_FFD2_3_3.txt')

    pathOut = 'design_out4.txt'
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
    run_optimization()