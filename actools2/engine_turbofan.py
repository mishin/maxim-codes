# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:07:43 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from flight_conditions import ISAtmosphere
from paths import MyPaths
from weight_tools import get_total_cg
from numpy import ones,array,linspace
from math import isnan
from engine_turbofan_analysis import engine_modeling
from scipy.interpolate import interp1d, Rbf
import matplotlib.pyplot as plt
import constants
from misc_tools import Normalization

class Propulsion(object):
    def __init__(self):
        self.numberOfEngines = 1
        self.engine = TurbofanEngine()
        self.CGx = None
        self.CGy = None
        self.CGz = None
        self.numberOfTanks = 1

    def load(self,name):
        self.engine.load(name)
        self._process_data()

    def _process_data(self):
        self.totalThrust = self.engine.thrustMC*self.numberOfEngines
        self._calc_cg()
        self._build_thrust_table()

    def _calc_cg(self):
        if hasattr(self.CGx,'__iter__'):
            mass = ones(self.numberOfEngines)*self.engine.mass
            self.totalMass, self.CG = get_total_cg(mass,self.CGx,self.CGy,self.CGz)


    def get_sfc(self,Mach,altitude,thrustReq):
        Mach      = float(Mach)
        altitude  = float(altitude)
        thrustReq = float(thrustReq) / self.numberOfEngines # thrust per engine
        
        sfcPerEngine = self.engine.get_sfc(Mach,altitude,thrustReq)
        return sfcPerEngine*self.numberOfEngines /3600. # kg/(N*sec)
    
    def get_sfc2(self,Mach,altitude,thrustReq):
        Mach = self._normMach.normalize(Mach)
        alt = self._normAlt.normalize(altitude)
        Treq = self._normTreq.normalize(thrustReq)
        return self.sfcModel(Mach,alt,Treq)
    
    def _build_thrust_table(self):
        n = 15
        MachList = linspace(0.05, 1.0, n)
        altList  = linspace(0, 2e4, n)
        Tmc = self.totalThrust
        TreqList = linspace(0.05*Tmc, Tmc, n)
        nsfc = n*n*n
        sfc       = ones(nsfc)
        Mach      = ones(nsfc)
        altitude  = ones(nsfc)
        thrustReq = ones(nsfc)
        
        self._normMach = Normalization(MachList[0], MachList[-1],0,1)
        self._normAlt  = Normalization(altList[0],altList[-1],0,1)
        self._normTreq = Normalization(TreqList[0],TreqList[-1],0,1)
        
        i = 0
        for M in MachList:
            for alt in altList:
                for Treq in TreqList:
                    _sfc = self.get_sfc(M,alt,Treq)
                    if not isnan(_sfc):
                        sfc[i] = _sfc
                        Mach[i] = self._normMach.normalize(M)
                        altitude[i] = self._normAlt.normalize(alt)
                        thrustReq[i] = self._normTreq.normalize(Treq)
                        i += 1
                        #print '%d\t%.4f\t%.0f\t%.0f\t%.6f'%(i, M, alt, Treq, _sfc)

        self.sfcModel = Rbf(Mach[:i], altitude[:i], thrustReq[:i], sfc[:i])

    def __call__(self,Mach,altitude,powerSetting):
        """ NOTE: powerSetting is linear function - will be replaced """
        thrustReq = self.engine.thrustMC * float(powerSetting)/100.
        return self.get_sfc(Mach,altitude,thrustReq)
    
    def get_thrust_available(self,altitude):
        atm = ISAtmosphere(altitude)
        rho0 = 1.2255
        return atm.density/rho0 *self.totalThrust

    def test(self):
        plt.figure()
        plt.hold(True)
        legend = list()
        alt = 1e4
        Mach = linspace(0.1,1,100)
        Tmc = self.totalThrust
        Pset = linspace(0.05*Tmc, Tmc ,5)
        print Pset
        sfc = ones([len(Pset),len(Mach)])
        sfc2 = ones([len(Pset),len(Mach)])

        for i,p in enumerate(Pset):
            for j,M in enumerate(Mach):
                sfc[i,j] = self.get_sfc(M,alt,p)
                sfc2[i,j] = self.get_sfc2(M,alt,p)

            plt.plot(Mach,sfc[i],marker='*')
            plt.plot(Mach,sfc2[i])
            
            legend.append('%.0f'%p)
            legend.append('%.0frbf'%p)
        plt.legend(legend)
        plt.show()


class TurbofanEngine(object):
    def __init__(self):
        self.name         = None
        self.thrustMC     = 0.0
        self.sfcMC        = 0.0
        self.mass         = 0.0
        self.length       = 0.0
        self.diameter     = 0.0
        # for thrust analysis
        self.designMach    = 0.0
        self.deignAltitude = 0.0
        self.designThrust  = 0.0
        T = array([362.873896304, 2267.9618519, 22679.618519])
        ratio = array([3.4,5.5,8])
        self._massCurve = interp1d(T,ratio,'linear')
        # for thrust table
        
    
    def load(self,name, xlsPath=None):
        """ loads turbofan engine from database file"""
        path = MyPaths()
        if xlsPath==None:
            xlsPath = path.db.engineTurbofan
        self.name = str(name)
        db = ReadDatabase(xlsPath,self.name)
        self.thrustMC       = db.read_row(-1,1) *constants.GRAVITY_ACCEL
        self.sfcMC          = db.read_row(-1,1)
        self.length         = db.read_row(-1,1)
        self.diameter       = db.read_row(-1,1)
        self.mass           = db.read_row(-1,1,False)
        self._calc_data()
    
    def _calc_data(self):
        if self.mass==0.0:
            self.mass = self._massCurve(self.thrustMC)
        if self.designThrust==0.0:
            self.designThrust = self.thrustMC
    
    def display(self):
        print self.name
        print self.thrustMC
        print self.mass
    
    def get_sfc(self,Mach,altitude,thrustReq):
        Tstatic = self.thrustMC
        Mcr = self.designMach
        altD = self.deignAltitude
        d, Td, sfcD, sfcF = engine_modeling(Tstatic,Mcr,altD,thrustReq,altitude,Mach)
        return sfcF


def run_test3():
    prop = Propulsion()
    prop.load('F404')
    prop.test()

if __name__=="__main__":
    run_test3()
        