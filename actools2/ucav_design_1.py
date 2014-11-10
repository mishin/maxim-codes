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
import sys

from misc_tools import read_tabulated_data_without_header, Normalization, RbfMod


class DesignFormulation(design.Design):
    def setup(self):
        self.neval = 0
        self.lb = np.array([40, 40, 6.00, 3.00, 0.5, 1.00, 3.000, -4, -4])
        #self.lb = np.array([40, 40, 6.00, 3.00, 0.5, 1.00, 3.000, -4, -4])
        self.ub = np.array([60, 60, 7.50, 5.25, 1.8, 1.80, 3.200,  0,  0])
        #self.ub = np.array([60, 60, 7.50, 5.25, 1.8, 1.80, 3.200,  0,  0])
        self.x0 = np.array([55, 55, 6.91, 4.15, 1.1, 1.44, 3.115,  0, -3])
        self.norm = Normalization(self.lb, self.ub,-1.,1.)
        self.xCurrent = np.zeros(len(self.x0))
        self.x0norm = self.norm.normalize(self.x0)
        #self.bnds = np.array([[l,u] for l,u in zip(self.lb,self.ub)])
        # --- constraints ---
        self.WemptyMax      = 3500.0
        self.CnbMin         = 0.003
        self.ClbMax         = -0.075 #-6.0e-4
        self.SMmin          = 0.05
        self.SMmax          = 0.15
        self.combatRadiusMin= 750.0
        self.RCmin          = 125.0
        self.VmaxMin        = 0.90 # Mach
        self.minSectionLen  = 1.4*self.propulsion.engine.length
        self.CmdeReq        = -0.01 # baseline -0.00866
        self.alphaTrimMax   = 8.0+1e-4 # deg
        self.elevTrimMax    = 20.0
        self.elevTrimMin    = -20.0
        self.gvfmAero = True
        self.Vtrim = 65.0 # 135kts: assumed value
        self._load_rbf_models()
        # --- payload ---
        self.SDB = MassComponent('drop payload', 1132.0, np.array([4.5, 0.0, 0.12]))
        # --- misc ---
        self.analysisData = np.zeros(10)
#        self._upd_approximation()
        # --- initial data ---
        self.Wf  = self.mass.fuel.mass
        self.CGf = self.mass.fuel.coords
        self.engineOffset = 0.75*self.propulsion.engine.length
        self.payloadOffset = 0.75*self.propulsion.engine.length
        self.fileFeasible = 'ucav_feasible.txt'
        self.set_x(self.x0norm)        
    
    def get_bounds(self,x0,c1=0.3):
        n = len(self.x0norm)
        #c1 = 0.25
        treshold = 0.005
        xUold = np.ones(n)
        xLold = -np.ones(n)
        xUnew = np.zeros(n)
        xLnew = np.zeros(n)
        u = xUold-x0
        l = x0 - xLold
        delta = np.zeros(len(self.x0))
        i=0
        for _l,_u in zip(l,u):
            delta[i] = min([_l,_u])*c1
            if delta[i]<treshold:
                d = c1*(xUold[i]-xLold[i])
                if _u<treshold:
                    xLnew[i] = x0[i] - d
                    xUnew[i] = x0[i]
                else:
                    xLnew[i] = x0[i]
                    xUnew[i] = x0[i] + d
            else:
                xUnew[i] = x0[i] + delta[i]
                xLnew[i] = x0[i] - delta[i]
            print '%.4f\t%.4f\t%.4f'%(_l,_u,delta[i])
            i += 1
        bounds = np.transpose(np.vstack([xLnew,xUnew]))
        return bounds
        
#        
#    def _upd_approximation(self):
#        pathIn = 'design_out4.txt'
#        pathInSamples = 'DOE_LHS250_FFD2_3_3.txt'
#        learnData = read_tabulated_data_without_header(pathIn)
#        xNorm = read_tabulated_data_without_header(pathInSamples)
#        learnData = np.transpose(learnData)
#        _Cnb     = RbfMod(xNorm, (learnData[2])*1e3)
#        self.Cnb = lambda x: _Cnb(x)/1e3
#        self.Clb = RbfMod(xNorm, learnData[3])
#        self.LD  = RbfMod(xNorm, learnData[0])

    def set_x(self,xnorm,fullUpdate=True):
        """
        To minimize function evaluation, this function re-calculates analysis if 
        new x is given, otherwise precalculated results are used.
        """
        self.xcurNorm = xnorm
        x = self.norm.denormalize(xnorm)
        try:
            self._upd_configuration(x)
            self._upd_analysis(xnorm,fullUpdate)
        except:
            print x
            print xnorm
            self.display_2d()
            self._upd_configuration(x)
            print 'Error occured'
            self._upd_analysis(xnorm,fullUpdate)
#            self.analysisData[0] = 0.0
#            self.analysisData[1:] = -np.ones(len(self.analysisData)-1)*100.0
#            raise ValueError
    
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
        cgNew2 = x+self.payloadOffset
        #self.mass.empty.update_item_cg('engine',cgNew,0,0)
        self.set_engine_cg(cgNew,0,0)
        self.SDB = MassComponent('drop payload', 1132.0, np.array([cgNew2, 0.0, 0.12]))
        #self.mass.payload.update_item_cg('drop payload',cgNew,.0,.0)
        
    
    def _upd_analysis(self,xnorm,fullUpdate):
        #self._update_mass()
        self._upd_drag()
        alt = self.designGoals.cruiseAltitude
        aero = self.get_aero_single_point(0.7,1e4,2.0)
        self.analysisData[0] = aero.LD

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
            self.analysisData[2] = aero.derivs.Cnb
            self.analysisData[3] = aero.derivs.Clb  #self.aero.derivs.Clb
            self.analysisData[4] = aero.SM
            self.analysisData[5] = self.combatRadius
            self.analysisData[7] = slf.run_max_TAS(alt).Mach
            aeroTrim = self.get_aero_trim(self.Vtrim, 0.0)
            self.analysisData[8] = aeroTrim.elevator
            self.analysisData[9] = aeroTrim.alpha

            self.mass.set_fuel_mass(self.Wf,self.CGf)
            if not self.mass.payload.item_exists(self.SDB.name):
                self.mass.payload.add_component(self.SDB)
            clm = ClimbDescent(self)
            self.analysisData[6] = clm.run_max_climb_rate(0).climbRate
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
        g[9] = self.analysisData[8] - self.elevTrimMin
        g[10] = _x[0] - _x[1]
        g[11] = self.alphaTrimMax - self.analysisData[9]
        print g>0
        out = ''
        for v in g:
            out += '%.0e '%v
        print out
        print '%.4f\t%.4f\t%.4f'%(self.analysisData[1],self.analysisData[5],self.analysisData[8])

        g[0] *= 1e-2
        g[1] *= 1e4
        g[2] *= 1e3
        g[3] *= 1e3
        g[4] *= 1e3
        g[5] *= 1e-1
        g[6] *= 1e-1
        g[7] *= 1e2
        g[9] *= 10.
        g[11] *= 10
        return g*1000.
    
    def get_aero_single_point(self,velocity,altitude,alpha=0.0,beta=0.0,
                          elevator=0.0,mass=None,cg=None,inertia=None,CD0=None):
        aero = super(DesignFormulation,self).get_aero_single_point(velocity,altitude,alpha,beta,elevator,mass,cg,inertia,CD0)
        if self.gvfmAero:
            sys.stdout.write('GVFM add> ')
            xcurrent = self.xcurNorm
            LDold = aero.coef.CL/(aero.coef.CD-self.CD0rbf(xcurrent))
            aero.coef.CL += self.CLrbf(xcurrent)
            aero.coef.CD += self.CDrbf(xcurrent)-self.CD0rbf(xcurrent)
            aero.LD      = LDold + self.LDrbf(xcurrent)
            aero.SM      += self.SMrbf(xcurrent)
            aero.CD0     += self.CD0rbf(xcurrent)
            aero.k       += self.krbf(xcurrent)
        return aero
    
    def update_aero_trim(self,velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                      mass=None,cg=None,inertia=None,CD0=None):
        super(DesignFormulation,self).update_aero_trim(velocity,altitude,CmTrim,loadFactor,mass,cg,inertia,CD0)
        if self.gvfmAero:
            xcurrent = self.xcurNorm
            #self.aeroResults.CD0     += self.CD0rbf(xcurrent)
            self.aeroResults.k       += self.krbf(xcurrent)

    def get_drag(self,velocity=None,altitude=None):
        cd0 = super(DesignFormulation,self).get_drag(velocity,altitude)
        if self.gvfmAero:
            return cd0 + self.CD0rbf(self.xcurNorm)
        else:
            return cd0

    def _load_rbf_models(self,DOEpath=None,rsltPath=None):

        if DOEpath==None:
            DOEpath = r'E:\1. Current work\2014 - UAV performance code\Results\DOE\LHS_dvar9_sample24.txt'
        if rsltPath==None:
            rsltPath = r'E:\1. Current work\2014 - UAV performance code\Results\DOE\CFDresults.txt'
        dCD0 = 0.0005
        xDOE = read_tabulated_data_without_header(DOEpath)
        fDOE = read_tabulated_data_without_header(rsltPath,1)
        self.CLrbf = RbfMod(xDOE, fDOE[:,4]-fDOE[:,13])
        self.CDrbf = RbfMod(xDOE, fDOE[:,5]-fDOE[:,14])
        self.CD0rbf= RbfMod(xDOE, fDOE[:,11]-fDOE[:,19]+dCD0)
        self.krbf  = RbfMod(xDOE, fDOE[:,12]-fDOE[:,20])
        self.SMrbf = RbfMod(xDOE, fDOE[:,10]-fDOE[:,18],0.5)
        self.LDrbf = RbfMod(xDOE, fDOE[:,7] -fDOE[:,15])
        self._xdoe = xDOE
    
    def show_results(self):
        print '{0:10}={1:10.4f}'.format('LD',self.analysisData[0])
        print '{0:10}={1:10.4f}'.format('SM',self.analysisData[4])
        print '{0:10}={1:10.4f}'.format('Cnb',self.analysisData[2])
        print '{0:10}={1:10.4f}'.format('Clb',self.analysisData[3])
        print '{0:10}={1:10.4f}'.format('elev',self.analysisData[8])
        print '{0:10}={1:10.4f}'.format('alpha',self.analysisData[9])
        print '{0:10}={1:10.4f}'.format('We',self.analysisData[1])
        print '{0:10}={1:10.4f}'.format('Rcombat',self.analysisData[5])
        print '{0:10}={1:10.4f}'.format('ROC',self.analysisData[6])
        print '{0:10}={1:10.4f}'.format('Mmax',self.analysisData[7])

def calculate_scaling_factors():
    def process_cd0(CL0=-.1,f=None):
        if f.ndim==2:
            CD1 = f[:,1]
            CD2 = f[:,4]
            CL1 = f[:,0]
            CL2 = f[:,3]
        else:
            CD1 = f[1]
            CD2 = f[4]
            CL1 = f[0]
            CL2 = f[3]
        tmp1 = (CL1-CL0)**2.
        tmp2 = (CL2-CL0)**2.
        k = (CD1 - CD2)/(tmp1 - tmp2)
        CD0 = CD1 - k*CL1*CL1
        return CD0, k
    
    def adjust_cd0(f,cd0avl,kavl):
        bnds = ((-.1,.1),)
        def obj(x,f):
            cd0,k = process_cd0(x,f)
            return 1e4*(cd0-cd0avl)**2.0 + (k-kavl)**2.0
        cl0 = fmin_slsqp(obj,[0],bounds=bnds,args=(f,))[0]
        cd0,k = process_cd0(cl0,f)
        return cd0,k,cl0

    xpath = r'E:\1. Current work\2014 - UAV performance code\Results\DOE\LHS_dvar9_sample24.txt'
    fpath = r'E:\1. Current work\2014 - UAV performance code\Results\DOE\CFD_FW_results-longitudinal2.txt'
    out = r'E:\1. Current work\2014 - UAV performance code\Results\DOE\CFD_FW_results-longitudinal3.txt'
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    xData = read_tabulated_data_without_header(xpath)
    fData = read_tabulated_data_without_header(fpath,1)
    
    k = np.zeros(fData.shape[0])
    CD0 = np.zeros(fData.shape[0])
    
    #CD0cfd, kcfd = process_cd0(f=fData)
    
    fid = open(out,'wt')
    fid.write('CD0avl\tkavl\tCD0cfd\tkcfd\tCD0add\tkadd\tCL0\n')
    for i,x in enumerate(xData):
        ac.set_x(x)
        aero = ac.get_aero_single_point(.7,1e4,alpha=2.)
        k[i] = aero.k
        CD0[i] = aero.CD0
        CD0cfd,kcfd,cl0 = adjust_cd0(fData[i],CD0[i],k[i])
        print CD0cfd, kcfd, cl0
        print 'case ',i
        fid.write('%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n'%(CD0[i],k[i],CD0cfd,kcfd,CD0cfd-CD0[i],kcfd-k[i],cl0))
    fid.close()
    
def run_optimization():
    print 'GVFM 1st iter: 24+1+1 points'
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    ac.set_x(ac.x0norm)
    for xval in ac.x0norm:
        sys.stdout.write('%.8f\t'%xval)
    sys.stdout.write('\n')

    #ac.display()
    
#    bnds = np.ones([len(ac.x0),2])
#    bnds[:,0] = -bnds[:,0]
    
    #low-fi
    #x0norm2 = np.array([ 0.20520747,-0.39226611,-0.47326701,-1.,0.14291925,0.98650447, 
    #                  0.37346922,  0.37756372,  0.65654722]) # low fidelity optimum
    x0norm2 = np.array([-0.10127195,-0.10127187,-0.60597138,-0.99999966,-0.71421434,0.97300896,1.00000803,0.57455958,0.3130853]) #gvfm 1
    #x0norm2 = np.array([-0.12782713,-0.12782713,-0.53160625,-0.99999966,-0.57132151,0.95951344,1.00000803,0.44481799,0.65654265])    
    x0 = ac._xdoe[25]
    bnds = ac.get_bounds(x0,0.25)
    print bnds

    rslt = fmin_slsqp(ac.f, x0, f_ieqcons=ac.g, bounds=bnds, iprint=2,epsilon=1e-2, acc=0.01)

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
    print 'alphaTrim',ac.analysisData[9]
    print 'elevTrim',ac.analysisData[8]
    #ac.aero.display()
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