# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:07:43 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from paths import MyPaths
from weight_tools import get_total_cg
from numpy import ones,array
from engine_turbofan_analysis import engine_modeling
from scipy.interpolate import interp1d


class Propulsion(object):
    def __init__(self):
        self.numberOfEngines = 0
        self.engine = TurbofanEngine()
        self.CGx = None
        self.CGy = None
        self.CGz = None
        self.numberOfTanks = 1
    def _process_data(self):
        self.totalThrust = self.engine.thrustMC*self.numberOfEngines
        self._calc_cg()
    def _calc_cg(self):
        if hasattr(self.CGx,'__iter__'):
            mass = ones(self.numberOfEngines)*self.engine.mass
            self.totalMass, self.CG = get_total_cg(mass,self.CGx,self.CGy,self.CGz)
    def get_sfc(self,Mach,altitude,thrustReq):
        Mach = float(Mach)
        altitude = float(altitude)
        thrustReq = float(thrustReq) / self.numberOfEngines # thrust per engine
        sfcPerEngine = self.engine.get_sfc(Mach,altitude,thrustReq)
        return sfcPerEngine*self.numberOfEngines
        
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
    
    def load(self,name, xlsPath=None):
        """ loads turbofan engine from database file"""
        path = MyPaths()
        if xlsPath==None:
            xlsPath = path.db.engineTurbofan
        self.name = str(name)
        db = ReadDatabase(xlsPath,self.name)
        self.thrustMC       = db.read_row(-1,1)
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
        Tstatic = self.designThrust
        Mcr = self.designMach
        altD = self.deignAltitude
        d, Td, sfcD, sfcF = engine_modeling(Tstatic,Mcr,altD,thrustReq,altitude,Mach)
        return sfcF


def run_test1():
    engine = TurbofanEngine()
    engine.load('F110')
    engine.display()

def run_test2():
    engine = TurbofanEngine()
    engine.load('newEngine1')
    engine.deignAltitude = 10000
    engine.designMach = 0.9
    engine.designThrust = engine.thrustMC
    
    alt = 8000.
    Mach = array([0.1,0.3,0.5,0.7,0.9])
    sfc = ones(len(Mach))
    for i,M in enumerate(Mach):
        sfc[i] = engine.get_sfc(M,alt,3500)
        print M,sfc[i]
    
if __name__=="__main__":
    run_test2()
        