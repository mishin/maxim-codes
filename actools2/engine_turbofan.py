# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:07:43 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from paths import MyPaths
from weight_tools import get_total_cg
from numpy import ones

class Propulsion(object):
    def __init__(self):
        self.numberOfEngines = 0
        self.engine = TurbofanEngine()
        self.CGx = None
        self.CGy = None
        self.CGz = None
        self.numberOfTanks = 1
    def _process_data(self):
        self.totalThrust = self.engine.thrustForced*self.numberOfEngines
        self._calc_cg()
    def _calc_cg(self):
        if hasattr(self.CGx,'__iter__'):
            mass = ones(self.numberOfEngines)*self.engine.mass
            self.totalMass, self.CG = get_total_cg(mass,self.CGx,self.CGy,self.CGz)

class TurbofanEngine(object):
    def __init__(self):
        self.name         = None
        self.thrustMC     = 0.0
        self.thrustForced = 0.0
        self.sfcMC        = 0.0
        self.scfForced    = 0.0
        self.mass         = 0.0
        self.length       = 0.0
        self.diameter     = 0.0
    
    def load(self,name, xlsPath=None):
        """ loads turbofan engine from database file"""
        path = MyPaths()
        if xlsPath==None:
            xlsPath = path.db.engineTurbofan
        self.name = str(name)
        db = ReadDatabase(xlsPath,self.name)
        self.thrustMC       = db.read_row(-1,1)
        self.thrustForced   = db.read_row(-1,1)
        self.sfcMC          = db.read_row(-1,1)
        self.scfForced      = db.read_row(-1,1)
        self.length         = db.read_row(-1,1)
        self.diameter       = db.read_row(-1,1)
        self.mass           = db.read_row(-1,1,False)
    
    def display(self):
        print self.name
        print self.thrustForced
        print self.mass
    
    def get_sfc(self,Mach,altitude,thrustReq):
        pass


def run_test1():
    engine = TurbofanEngine()
    engine.load('F110')
    engine.display()

if __name__=="__main__":
    run_test1()
        