# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 21:43:04 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from paths import MyPaths
import numpy as np
import airfoil
from engine_turbofan import Propulsion
from weight_fw import get_flying_wing_mass, AircraftMass
from flight_conditions import FlightConditions
from display_aircraft import flying_wing_display
from drag import get_friction_drag_FW

path = MyPaths()

class FlyingWing(object):
    def __init__(self):
        self.wing = Wing() #main geometry
        self.designGoals = DesignGoals()
        self.landingGear = LandingGear()
        self.vlm = VLMparameters()
        self.weight = AircraftMass() #total mass list
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
        self.weight.fuel.add_item('Fuel tank Right', self.designGoals.fuelMass/2.0)
        self.weight.fuel.add_item('Fuel tank Left', self.designGoals.fuelMass/2.0)
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
            self.weight.payload.add_item(name,mass,CG,MOI)
        self.weight.update_total()
        self._process_data()
        
        #TODO: tmp assumptions ---
        self.propulsion.totalThrust = 13000.
        self.propulsion.engineMass = 1000.
        self.designGoals._process_data()

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

    def _update_weight(self):
        self.weight = get_flying_wing_mass(self)
    
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

class Wing(object):
    def __init__(self):
        self.locationLE        = np.zeros(3)
        self.segSpans          = None
        self.chords            = None
        self.airfoilNames      = None
        self.airfoils          = list()
        self.secOffset         = None
        self.segDihedral       = None
        self.segAreas          = None
        self.secTwist          = None
        self.incidence         = 0.0
        self.secAngles         = None
        self.leadingEdge       = np.zeros(3)
        self.csAileronLocation = None
        self.csFlapLocation    = None
        self.csFlapType        = None
        self.material          = None
        self.MAC               = 0.0
        self.MAClocation       = np.zeros(2) # x,y
        self.nSeg              = 0
        self.secThickness      = None

    def _process_data(self):
        self.nSeg = len(self.segSpans)
        self.nSec = self.nSeg+1
        self._load_airfoils()
        self._calc_geometry_data()

    def _load_airfoils(self,xlsPath=None):
        if xlsPath==None:
            xlsPath = path.db.airfoil
        self.secThickness = np.zeros(self.nSec)
        for i,name in enumerate(self.airfoilNames):
            self.airfoils.append(airfoil.load(name))
            self.secThickness[i] = self.airfoils[i].thickness

    def _calc_geometry_data(self):
        self._calc_apex()
        self._calc_mac()
        self._calc_span()
        self._calc_wetted_area()
        self._calc_angles()
        self._calc_elastic_axis_sweep()
        self._calc_controls()
    
    def _calc_controls(self):
        self.aileronArea = self._get_cs_area(self.csAileronLocation)
        self.flapArea = self._get_cs_area(self.csFlapLocation)
        self.csArea = self.aileronArea + self.flapArea
    
    def _get_cs_area(self,location):
        startSec = np.argmax(location<1)
        endSec = len(location) - np.argmax(location[::-1]<1)-1
        area = 0.0
        for i in range(startSec,endSec):
            c1 = self.chords[i]*(1-location[i])
            c2 = self.chords[i+1]*(1-location[i+1])
            area += (c1 + c2)*self.segSpans[i]
        return area

    def _calc_span(self):
        self.span = 2.0*np.sum(self.segSpans)
        self.aspectRatio = self.span**2.0/self.area

    def _calc_apex(self):
        # calculate leading edge point of each section - section incidence is 
        # not considered
        secApex = np.zeros([3,self.nSec])
        secApex[0,1:] = self.secOffset.cumsum()
        secApex[1,1:] = self.segSpans.cumsum()
        secApex[2,1:] = secApex[0,1:]*np.tan(np.radians(self.segDihedral))
        self.secApex = np.transpose(secApex)
    
    def _calc_angles(self):
        # calculates angle of each section
        self.secAngles = np.zeros(self.nSec)
        self.secAngles[0] = self.incidence
        for i,twist in enumerate(self.secTwist):
            self.secAngles[i+1] = self.secAngles[i]+twist
    
    def _calc_elastic_axis_sweep(self,FSfrontSpar=0.3,FSaftSpar=0.7):
        """
        Assumed values are location of front and aft spars. Typically the elastic 
        axis is located between two spars. 
        NOTE: This calculation is very approximate.
        """
        eaLoc = (FSaftSpar-FSfrontSpar)/2.0
        xr = self.secApex[0,0] + self.chords[0]*eaLoc
        xt = self.secApex[-1,0] + self.chords[-1]*eaLoc
        yr = self.secApex[0,1]
        yt = self.secApex[-1,1]
        self.sweepElasticRad = np.arctan((yt-yr)/(xt-xr))
        self.sweepElasticDeg = np.degrees(self.sweepElasticRad)
        dx, dy = xt-xr, yt-yr
        self.structuralSpan = (dx*dx+dy*dy)**0.5

    def _calc_mac(self):
        """ calculates mean aerodynamic chord and it's location"""
        self.segMAC = np.zeros(self.nSeg)
        self.segAreas = np.zeros(self.nSeg)
        def mac_of_one_segment(c1, c2, x1, x2, y1, y2):
            span = y2 - y1
            area = 0.5*(c1+c2)*span
            taper = c2/c1
            taper1 = taper + 1.0
            frac = (taper+taper1)/(3.0*taper1)
            mac = c1*(taper*taper + taper + 1.0) / (1.5*taper1)
            ymac = y1 + frac*span
            xlemac = x1 + frac*(x2-x1)
            return area, mac, ymac, xlemac
        area, mac, ymac, xlemac = 0.0, 0.0, 0.0, 0.0
        for i in range(self.nSeg):
            c1,c2 = self.chords[i], self.chords[i+1]
            x1,x2 = self.secApex[i,0], self.secApex[i+1,0]
            y1,y2 = self.secApex[i,1], self.secApex[i+1,1]
            aseg, macseg, yseg, xseg = mac_of_one_segment(c1,c2,x1,x2,y1,y2)
            self.segMAC[i] = macseg
            self.segAreas[i] = aseg
            area   += aseg
            mac    += macseg*aseg
            ymac   += yseg*aseg
            xlemac += xseg*aseg
        self.MAC = mac/area
        self.MAClocation = np.array([xlemac/area, ymac/area])
        self.area = 2.0*area

    def _calc_wetted_area(self):
        wettedArea = 0.0
        for i in range(len(self.segSpans)-1):
            root = self.airfoils[i].length*self.chords[i]
            tip = self.airfoils[i+1].length*self.chords[i+1]
            wettedArea += (root+tip)*self.segSpans[i+1]
        self.wettedArea = wettedArea


class VLMparameters(object):
    def __init__(self):
        self.panelsChordwise = 0
        self.panelsSpanwise = 0
        self.distribution = None

def run_test2():
    ac = FlyingWing()
    ac.load_xls('sample_B45c')
    ac._update_weight()
    #ac.display()

def run_test1():
    import matplotlib.pyplot as plt
    ac = FlyingWing()
    ac.load_xls('sample1')
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