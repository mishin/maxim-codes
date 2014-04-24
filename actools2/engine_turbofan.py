# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:07:43 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from paths import MyPaths


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
        self.totalMass = self.engine.mass*self.numberOfEngines

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
        self.mass           = db.read_row(-1,1)
    
    def display(self):
        print self.name
        print self.thrustForced
        print self.mass


def run_test1():
    engine = TurbofanEngine()
    engine.load('F110')
    engine.display()

if __name__=="__main__":
    run_test1()
        