# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 21:43:04 2014

@author: Maxim
"""
from db_tools import ReadDatabase, WriteDatabase
from paths import MyPaths
import numpy as np
import airfoil
from engine_turbofan import TurbofanEngine
from weight import AircraftMass
from flight_conditions import FlightConditions

path = MyPaths()

class FlyingWing(object):
    def __init__(self):
        self.wing = Wing() #main geometry
        self.designGoals = DesignGoals()
        self.landingGear = LandingGear()
        self.vlm = VLMparameters()
        self.weight = AircraftMass() #total mass list
        self.drag = None #total drag list
        self.propulsion = TurbofanEngine() #table lookup database
        
    def load_xls(self, name, xlsPath=None):
        if xlsPath==None:
            xlsPath = path.db.aircraftFW
        keyword = 'SECTION: '
        db = ReadDatabase(xlsPath, name, keyword)
        idx = db.find_header(keyword+'DESIGN QUANTITIES')
        self.type                          = db.read_row(idx+1,1)
        self.designGoals.fuelMass          = db.read_row(-1,1)
        self.designGoals.grossMass         = db.read_row(-1,1)
        self.designGoals.cruiseSpeed       = db.read_row(-1,1)
        self.designGoals.CruiseAltitude    = db.read_row(-1,1)
        self.designGoals.loadFactor        = db.read_row(-1,1)
        self.designGoals.loadFactorLanding = db.read_row(-1,1)
        self.designGoals.numberOfOccupants = db.read_row(-1,1)
        self.weight.fuel.add_item('Fuel tank Right', self.designGoals.fuelMass/2.0)
        self.weight.fuel.add_item('Fuel tank Left', self.designGoals.fuelMass/2.0)
        idx = db.find_header(keyword+'MAIN WING')
        self.wing.segments           = db.read_row(idx+1,1,True)
        self.wing.chords             = db.read_row(-1,1,True)
        self.wing.airfoilNames       = db.read_row(-1,1,True)
        self.wing.secOffset          = db.read_row(-1,1,True)
        self.wing.segDihedral        = db.read_row(-1,1,True)
        self.wing.secTwist           = db.read_row(-1,1,True)
        self.wing.incidence          = db.read_row(-1,1,False)
        self.wing.csAileronLocation  = db.read_row(-1,1,True)
        self.wing.csFlapLocation     = db.read_row(-1,1,True)
        self.wing.csFlapType         = db.read_row(-1,1,False)
        self.wing.material           = db.read_row(-1,1,False)
        idx = db.find_header(keyword+'PROPULSION')
        self.propulsion.numberOfEngines = db.read_row(idx+1,1,True)
        self.propulsion.CGx             = db.read_row(-1,1,True)
        self.propulsion.CGy             = db.read_row(-1,1,True)
        self.propulsion.CGz             = db.read_row(-1,1,True)
        self.propulsion.mass            = db.read_row(-1,1)
        self.propulsion.tableName       = db.read_row(-1,1)
        idx = db.find_header(keyword+'LANDING GEAR')
        self.landingGear.groundContactX = db.read_row(idx+1,1,True)
        self.landingGear.groundContactY = db.read_row(-1,1,True)
        self.landingGear.groundContactZ = db.read_row(-1,1,True)
        self.landingGear.tireWidth      = db.read_row(-1,1,True)
        self.landingGear.tireDiameter   = db.read_row(-1,1,True)
        idx = db.find_header(keyword+'VLM PARAM')
        self.vlm.panelsChordwise = db.read_row(idx+1,1,False)
        self.vlm.panelsSpanwise  = db.read_row(-1,1,False)
        self.vlm.distribution    = db.read_row(-1,1,False)
        sec = db.read_section('PAYLOAD',0)
        for line in sec:
            name = str(line[0])
            mass = float(line[1])
            CG   = np.array([float(val) for val in line[2:5]])
            MOI  = np.array([float(val) for val in line[5:8]])
            self.weight.payload.add_item(name,mass,CG,MOI)
        self.weight.update_total()
        self._process_data()
    
    def _process_data(self):
        self.wing._process_data()
        self.designGoals.set_flight_conditions(self.wing.MAC)
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
    def set_flight_conditions(self, refLength=1.0):
        self.fc = FlightConditions(self.cruiseSpeed, self.cruiseAltitude, refLength)

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
        self.MAC               = 0.0
        self.MAClocation       = np.array([0,0,0])

    def _process_data(self):
        self._load_airfoils()
        self._calc_geometry_data()

    def _load_airfoils(self,xlsPath=None):
        if xlsPath==None:
            xlsPath = path.db.airfoil
        for name in self.airfoilNames:
            self.airfoils.append(airfoil.load(name))

    def _calc_geometry_data(self):
        self._calc_mac()
        self._calc_projected_area()
        self._calc_wetted_area()

    def _calc_mac(self):
        for i,segSpan in enumerate(self.segments[1:]):
            pass

    def _calc_wetted_area(self):
        n = len(self.segments)-1
        area = np.zeros(n)
        for i in range(n):
            root = self.airfoils[i].length*self.chords[i]
            tip = self.airfoils[i+1].length*self.chords[i+1]
            area[i] = (root+tip)*self.segments[i+1]
        self.wettedArea = area.sum()
    def _calc_projected_area(self):
        segArea = self.segments[1:]*(self.chords[:-1]+self.chords[1:])
        self.area = segArea.sum()

class VLMparameters(object):
    def __init__(self):
        self.panelsChordwise = 0
        self.panelsSpanwise = 0
        self.distribution = None

def run_test1():
    ac = FlyingWing()
    ac.load_xls('sample1')
    print ac.wing.wettedArea

if __name__=="__main__":
    run_test1()