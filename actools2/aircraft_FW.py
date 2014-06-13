# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 21:43:04 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from paths import MyPaths
import numpy as np
from scipy.interpolate import Akima1DInterpolator

from engine_turbofan import Propulsion
from weight_fw import get_flying_wing_mass, AircraftMass
from flight_conditions import FlightConditions
from display_aircraft import flying_wing_display
from drag import get_friction_drag_FW
from wing import Wing
from aero_avl_fw import Aerodynamics, FlightConditionsAVL

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
        self.aeroResults = None
        
        # FIXME: drag calculation is temporal - DO NOT USE for real calculation
        cd = np.array([0.01718,0.01718, 0.01662, 0.01622, 0.01589, 0.0156, 0.01532, 0.01506, 
                       0.0148, 0.03573, 0.04254, 0.03174, 0.02919, 0.02753, 0.02602, 
                       0.02469, 0.0235, 0.02245, 0.02149])
        M = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 
                      1.4, 1.5, 1.6, 1.7, 1.8, 1.9])
        self._dragCurve = Akima1DInterpolator(M,cd)

        
    def load_xls(self, name, xlsPath=None):
        self.name = str(name)
        if xlsPath==None:
            xlsPath = path.db.aircraftFW
        keyword = 'SECTION: '
        db = ReadDatabase(xlsPath, name, keyword)
        idx = db.find_header(keyword+'DESIGN QUANTITIES')
        self.type                          = db.read_row(idx+1,1)
        self.designGoals.fuelMass          = db.read_row(-1,1)
        self.designGoals.avionicsMass      = db.read_row(-1,1)
        self.designGoals.grossMass         = db.read_row(-1,1)
        self.designGoals.cruiseMach        = db.read_row(-1,1)
        self.designGoals.cruiseAltitude    = db.read_row(-1,1)
        self.designGoals.loadFactor        = db.read_row(-1,1)
        self.designGoals.loadFactorLanding = db.read_row(-1,1)
        self.designGoals.staticThrust      = db.read_row(-1,1)
        self.designGoals.numberOfOccupants = db.read_row(-1,1)
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
        elevonLocation               = db.read_row(-1,1,True)
        flapLocation                 = db.read_row(-1,1,True)
        self.wing.set_elevon(elevonLocation)
        self.wing.set_flap(flapLocation)
        self.wing.flap.type         = db.read_row(-1,1,False)
        self.wing.material           = db.read_row(-1,1,False)
        self.wing.fuelTankCGratio    = db.read_row(-1,1,True)
        idx = db.find_header(keyword+'PROPULSION')
        engineName                      = db.read_row(idx+1,1,False)
        self.propulsion.engine.load(engineName)
        self.propulsion.numberOfEngines = int(db.read_row(-1,1,False))
        self.propulsion.CGx             = db.read_row(-1,1,True)
        self.propulsion.CGy             = db.read_row(-1,1,True)
        self.propulsion.CGz             = db.read_row(-1,1,True)
        self.propulsion.engine.deignAltitude = self.designGoals.cruiseAltitude
        self.propulsion.engine.designMach    = self.designGoals.cruiseMach
        self.propulsion.engine.designThrust  = self.designGoals.staticThrust / self.propulsion.numberOfEngines
        idx = db.find_header(keyword+'LANDING GEAR')
        self.landingGear.groundContactX = db.read_row(idx+1,1,True)
        self.landingGear.groundContactY = db.read_row(-1,1,True)
        self.landingGear.groundContactZ = db.read_row(-1,1,True)
        self.landingGear.tireWidth      = db.read_row(-1,1,True)
        self.landingGear.tireDiameter   = db.read_row(-1,1,True)
        self.landingGear.strutLength    = db.read_row(-1,1,True)
        #self.landingGear._set_right_mlg()
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
    
    def update_aero_trim(self,velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                      mass=None,cg=None,inertia=None,CD0=None):
        aero = Aerodynamics(self)
        fc = FlightConditionsAVL(self,velocity,altitude,CmTrim,loadFactor,mass,
                                 cg,inertia,CD0)
        self.aeroResults = aero.run_trim(fc)
    
    def get_aero_trim(self,velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                      mass=None,cg=None,inertia=None,CD0=None):
        """ run avl analysis and return trim state """
        self.update_aero_trim(velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                              mass=None,cg=None,inertia=None,CD0=None)
        return self.aeroResults

    def get_cg(self,update=True):
        # SI units [x,y,z]
        if update:
            self._update_mass()
        return self.mass.total.CG

    def get_mass(self,update=True):
        if update:
            self._update_mass()
        return self.mass.total.totalMass

    def get_mass_empty(self,update=True):
        if update:
            self._update_mass()
        return self.mass.empty.totalMass

    def get_drag(self,velocity=None,altitude=None):
        if velocity==None:
            velocity = self.designGoals.cruiseSpeed
        if altitude==None:
            altitude = self.designGoals.cruiseAltitude
        fc = FlightConditions(velocity,altitude)
        cd = self._dragCurve(fc.Mach)
        return cd #FIXME: for debug only
    
    def get_inertia(self):
        return np.zeros(3) #FIXME: for debug only
    
    def get_sfc(self,velocity,altitude,thrustRequired):
        fc = FlightConditions(velocity,altitude,0.0,self.wing.MAC)
        sfc = self.propulsion.get_sfc(fc.Mach,altitude,thrustRequired)
        return sfc


class LandingGear(object):
    def __init__(self):
        self.tireWidth = None
        self.tireDiameter = None
        self.groundContactX = None
        self.groundContactY = None
        self.groundContactZ = None
        self.strutLength    = None

#    def _set_right_mlg(self):
#        assert len(self.groundContactX)>2
#        self.groundContactX = np.hstack([self.groundContactX,self.groundContactX[1]])
#        self.groundContactY = np.hstack([self.groundContactY,-self.groundContactY[1]])
#        self.groundContactZ = np.hstack([self.groundContactZ,self.groundContactZ[1]])
#        self.strutLength = np.hstack([self.strutLength,self.strutLength[1]])
#        self.tireWidth = np.hstack([self.tireWidth,self.tireWidth[1]])
#        self.tireDiameter = np.hstack([self.tireDiameter,self.tireDiameter[1]])
#    
    def get_cg(self):
        pass
    
    def get_cg_retracted(self):
        pass

class DesignGoals(object):
    def __init__(self):
        self.fuelMass          = 0.0
        self.avionicsMass      = 0.0
        self.grossMass         = 0.0
        self.cruiseSpeed       = 0.0
        self.cruiseAltitude    = 0.0
        self.loadFactor        = 0.0
        self.loadFactorLanding = 0.0
        self.staticThrust      = 0.0
        self.numberOfOccupants = 0
    def set_flight_conditions(self, refLength=1.0):
        self.fc = FlightConditions(self.cruiseMach, self.cruiseAltitude, refLength)
    def _process_data(self,MAC=1.0):
        self.set_flight_conditions(MAC)
        self.cruiseMach = self.fc.Mach
        self.cruiseSpeed = self.fc.velocity
    def __repr__(self):
        out = 'Design goals\n====================\n'
        for attr, value in self.__dict__.iteritems():
            if type(value) is float or type(value) is np.float64:
                out += '{:<15} = {:<12.2f}\n'.format(attr,value)
        return out


class VLMparameters(object):
    def __init__(self):
        self.panelsChordwise = 0
        self.panelsSpanwise  = 0
        self.distribution    = None

def run_test3():
    import matplotlib.pyplot as plt
    ac = FlyingWing()
    ac.load_xls('X47A')
    M = np.linspace(0,1.2,50)
    cd = np.array([ac.get_drag(m,1e4) for m in M])
    plt.plot(M,cd,'bs-')
    plt.show()

def run_test2():
    ac = FlyingWing()
    ac.load_xls('X47A')
    print ac.wing.area
    ac.mass.display()
    ac.get_aero_trim(150,3000)
    ac.aeroResults.display()

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

def run_test4():
    import matplotlib.pyplot as plt
    ac = load('X45C')

if __name__=="__main__":
    run_test3()