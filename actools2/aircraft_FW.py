# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 21:43:04 2014

@author: Maxim
"""
from db_tools import ReadDatabase, WriteDatabase
from paths import MyPaths
import numpy as np
import airfoil

path = MyPaths()

class FlyingWing(object):
    def __init__(self):
        self.wing = Wing() #main geometry
        self.designGoals = DesignGoals()
        self.landingGear = LandingGear()
        self.vlm = VLMparameters()
        self.mass = None #total mass list
        self.drag = None #total drag list
        self.propulsion = None #table lookup database
        
    def load_xls(self, name, xlsPath=None):
        if xlsPath==None:
            xlsPath = path.db.aircraftFW
        db = ReadDatabase(xlsPath, name, 'SECTION: ')
        sec = db.read_section('DESIGN QUANTITIES')
        self.type = sec[0]
        self.designGoals.fuelMass          = sec[1]
        self.designGoals.grossMass         = sec[2]
        self.designGoals.cruiseSpeed       = sec[3]
        self.designGoals.CruiseAltitude    = sec[4]
        self.designGoals.loadFactor        = sec[5]
        self.designGoals.loadFactorLanding = sec[6]
        self.designGoals.numberOfOccupants = sec[7]
        sec = db.read_section('MAIN WING')
        self.wing.segments           = sec[0]
        self.wing.chords             = sec[1]
        self.wing.airfoilNames       = sec[2]
        self.wing.secOffset          = sec[3]
        self.wing.segDihedral        = sec[4]
        self.wing.secTwist           = sec[5]
        self.wing.incidence          = sec[6]
        self.wing.csAileronLocation  = sec[7]
        self.wing.csFlapLocation     = sec[8]
        self.wing.csFlapType         = sec[9]
        self.wing.material           = sec[10]
        sec = db.read_section('LANDING GEAR')
        self.landingGear.groundContactX = sec[0]
        self.landingGear.groundContactY = sec[1]
        self.landingGear.groundContactZ = sec[2]
        self.landingGear.tireWidth      = sec[3]
        self.landingGear.tireDiameter   = sec[4]
        sec = db.read_section('VLM PARAM')
        self.vlm.panelsChordwise = sec[0]
        self.vlm.panelsSpanwise  = sec[1]
        self.vlm.distribution    = sec[2]
    def save_xls(self):
        pass


class LandingGear(object):
    def __init__(self):
        self.tireWidth = 0.0
        self.tireDiameter = 0.0
        self.groundContactX = None
        self.groundContactY = None
        self.groundContactZ = None

class DesignGoals(object):
    def __init__(self):
        self.fuelMass = 0.0
        self.grossMass = 0.0
        self.cruiseSpeed = 0.0
        self.cruiseAltitude = 0.0
        self.loadFactor = 0.0
        self.loadFactorLanding = 0.0
        self.numberOfOccupants = 0

class Wing(object):
    def __init__(self):
        self.segments          = None
        self.chords            = None
        self.airfoilNames      = None
        self.airfoils          = list()
        self.secOffset         = None
        self.segDihedral       = None
        self.secTwist          = None
        self.incidence         = 0.0
        self.leadingEdge       = np.zeros(3)
        self.csAileronLocation = None
        self.csFlapLocation    = None
        self.csFlapType        = None
        self.material          = None
    
    def load_airfoils(self,xlsPath=None):
        if xlsPath==None:
            xlsPath = path.db.airfoil
        for name in self.airfoilNames:
            self.airfoils.append(airfoil.load(name))
    def get_geometry_data(self):
        pass
        
class VLMparameters(object):
    def __init__(self):
        self.panelsChordwise = 0
        self.panelsSpanwise = 0
        self.distribution = None

def run_test1():
    ac = FlyingWing()
    ac.load_xls('sample1')
    print ac.designGoals.cruiseAltitude
    print ac.wing.chords
    print ac.wing.airfoilNames
    print ac.wing.material
    print ac.wing.csFlapLocation


if __name__=="__main__":
    run_test1()