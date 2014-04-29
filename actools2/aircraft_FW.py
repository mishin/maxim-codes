# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 21:43:04 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from paths import MyPaths
import numpy as np

from engine_turbofan import Propulsion
from weight_fw import get_flying_wing_mass, AircraftMass
from flight_conditions import FlightConditions
from display_aircraft import flying_wing_display
from drag import get_friction_drag_FW
from wing import Wing

path = MyPaths()

def load(name):
    ac = FlyingWing()
    ac.load_xls(name)
    return ac


class FlyingWing(object):
    def __init__(self):
        self.wing = Wing() #main geometry
        self.designGoals = DesignGoals()
        self.landingGear = LandingGear()
        self.vlm = VLMparameters()
        self.mass = AircraftMass() #total mass list
        self.drag = None #total drag list
        self.propulsion = Propulsion() #table lookup database
        
    def load_xls(self, name, xlsPath=None):
        self.name = str(name)
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
#        self.mass.fuel.add_item('Fuel tank Right', self.designGoals.fuelMass/2.0)
#        self.mass.fuel.add_item('Fuel tank Left', self.designGoals.fuelMass/2.0)
        idx = db.find_header(keyword+'FUSELAGE')
        self.fusWidth = db.read_row(idx+1,1,False)
        idx = db.find_header(keyword+'MAIN WING')
        self.wing.segSpans           = db.read_row(idx+1,1,True)
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
        self.wing.fuelTankCGratio    = db.read_row(-1,1,True)
        idx = db.find_header(keyword+'PROPULSION')
        engineName                      = db.read_row(idx+1,1,False)
        self.propulsion.engine.load(engineName)
        self.propulsion.numberOfEngines = db.read_row(-1,1,True)
        self.propulsion.CGx             = db.read_row(-1,1,True)
        self.propulsion.CGy             = db.read_row(-1,1,True)
        self.propulsion.CGz             = db.read_row(-1,1,True)
        idx = db.find_header(keyword+'LANDING GEAR')
        self.landingGear.groundContactX = db.read_row(idx+1,1,True)
        self.landingGear.groundContactY = db.read_row(-1,1,True)
        self.landingGear.groundContactZ = db.read_row(-1,1,True)
        self.landingGear.tireWidth      = db.read_row(-1,1,True)
        self.landingGear.tireDiameter   = db.read_row(-1,1,True)
        self.landingGear.strutLength    = db.read_row(-1,1,True)
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
            self.mass.payload.add_item(name,mass,CG,MOI)
        self.mass.update_total()
        self._process_data()
        fuelCG = self.wing.locate_on_wing(self.wing.fuelTankCGratio[0],self.wing.fuelTankCGratio[1])
        fuelCG2 = np.array([fuelCG[0],-fuelCG[1],fuelCG[2]])
        self.mass.fuel.add_item('Fuel tank right',self.designGoals.fuelMass/2.,fuelCG)
        self.mass.fuel.add_item('Fuel tank left',self.designGoals.fuelMass/2.,fuelCG2)
        #TODO: tmp assumptions ---
        self.propulsion.totalThrust = 13000.
        self.propulsion.engineMass = 1000.
        self.designGoals._process_data()
        self._update_mass()

    def _process_data(self):
        self.wing._process_data()
        self.designGoals.set_flight_conditions(self.wing.MAC)
        self.designGoals._process_data()
        self.propulsion._process_data()
    def save_xls(self):
        pass

    def display(self,showAxes=False):
        flying_wing_display(self,showAxes)
    
    def _update_parasite_drag(self,velocity,altitude):
        if velocity==None:
            velocity = self.designGoals.cruiseSpeed
        if altitude==None:
            altitude = self.designGoals.cruiseAltitude
        self.drag = get_friction_drag_FW(self,velocity,altitude)

    def _update_mass(self):
        self.mass = get_flying_wing_mass(self)
    
    def get_parasite_drag(self,velocity=None,altitude=None):
        """
        Returns friction and form drag coefficient at given velocity and altitude. 
        If velocity and altitude are not specified then values from design targets
        are used
        
        Parameters
        ----------
        
        velocity : float, m/sec
        altitude : float, m
        """
        self._update_parasite_drag(velocity,altitude)
        return self.drag.get_total_drag()

    def get_cg(self):
        pass
    def get_mass_total(self):
        pass
    def get_mass_empty(self):
        pass


class LandingGear(object):
    def __init__(self):
        self.tireWidth = 0.0
        self.tireDiameter = 0.0
        self.groundContactX = None
        self.groundContactY = None
        self.groundContactZ = None
        self.strutLength    = None

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
    def _process_data(self):
        self.cruiseMach = self.cruiseSpeed/self.fc.atm.soundSpeed

        






class VLMparameters(object):
    def __init__(self):
        self.panelsChordwise = 0
        self.panelsSpanwise = 0
        self.distribution = None

def run_test2():
    ac = FlyingWing()
    ac.load_xls('sample_B45c')
    #ac.mass.display()
    #ac.mass.airframe.display()
    #ac.display()
    print ac.wing.area
    print ac.wing.MAC
    print ac.wing.span
    print ac.wing.segSweepC2deg
    print ac.wing.equivSweepLEdeg
    print ac.wing.equivThickness
    print ac.wing.equivSweepC2deg
    print ac.wing.equivSweepC4deg

def run_test1():
    import matplotlib.pyplot as plt
    ac = FlyingWing()
    ac.load_xls('sample_B45c')
#    print ac.wing.MAC
#    print ac.wing.MAClocation
#    print ac.wing.secThickness
    
    h = ac.designGoals.cruiseAltitude
    Vcrs = ac.designGoals.cruiseSpeed
    V = np.linspace(0.1*Vcrs,Vcrs,20)
    a = ac.designGoals.fc.atm.soundSpeed
    rho = ac.designGoals.fc.atm.density
    cd = np.zeros(len(V))
    D = np.zeros(len(V))
    plt.figure(1)
    for i,v in enumerate(V):
        cd[i] = ac.get_parasite_drag(v,h)
        qS = rho*v*v/2.*ac.wing.area
        D[i] = qS*cd[i]
    plt.plot(V/a,cd,'bo-')
    plt.xlabel('Mach number')
    plt.ylabel('Friction+Profile drag coefficient')
    plt.grid(True)
    plt.show()
    ac.drag.display()
    ac.display()

if __name__=="__main__":
    run_test2()