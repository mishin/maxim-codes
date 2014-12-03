# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 21:43:04 2014
@author: Maxim

This is the main module that stores information about the aircraft and provides 
easy interface to the basic analysis modules through set of methods.

Examples
--------

Import module

>>> import aircraft_FW

Load aircraft configuration from Excel database and store it to variable *ac*

>>> ac = aircraft_FW.load('X45C')

Display 2D planform shape of an aircraft

>>> ac.display_2d()

Display aircraft in 3D using Mayavi

>>> ac.display()

Weight and balance
~~~~~~~~~~~~~~~~~~

Component mass and CG are calculated automatically when aircraft is loaded. 
Currently equations for fighter aircraft are used. Weight estimation method should 
be changed for high aspect ratio wings. Required values can be accessed through 
special methods of aircraft class as shown below:

Get total mass of an aircraft

>>> ac.get_mass()
7657.0258133023435

Get empty mass of an aircraft

>>> ac.get_mass_empty()
3625.025813302344

Get CG location

>>> ac.get_cg()
array([ 4.09514882,  0.        ,  0.03761437])

Object 'ac.mass' gives extended access to weight and balance parameters. Total 
mass, empty mass and CG can be also accessed through it:

>>> ac.mass.total()
7657.0258133023435

>>> ac.mass.empty()
3625.0258133

>>> ac.mass.total.CG
array([ 4.09514882,  0.        ,  0.03761437])

Mass and CG of airframe (empty) and payload can be accessed separately.

>>> ac.mass.empty.CG
array([ 4.12844815,  0.        ,  0.04983231])

>>> ac.mass.payload() # total payload mass
1032.0

>>> ac.mass.payload.CG
array([ 4.115,  0.   ,  0.12 ])

Display list of components

>>> ac.mass.empty.display()
===============================================================
X45C mass components breakdown
===============================================================
Name                |Mass,kg |CG:X,m    |CG:Y,m    |CG:Z,m    |
---------------------------------------------------------------
wing                 701.69     4.1333     0.0000     0.0000
fuselage             114.06     3.4550     0.0000     0.0000
main landing gear    277.61     5.0000     0.0000     0.0000
nose landing gear     48.99     1.0000     0.0000     0.0000
engine              1036.00     4.1150     0.0000     0.1000
engine mount          65.29     4.1150     0.0000     0.1000
engine section        15.67     4.1150     0.0000     0.1000
engine cooling        11.59     4.1150     0.0000     0.1000
air intake            95.69        N/A        N/A        N/A
oil cooling           17.15        N/A        N/A        N/A
engine controls        6.17     4.1150     0.0000     0.1000
engine starter        75.75        N/A        N/A        N/A
fuel system          253.23        N/A        N/A        N/A
flight controls       96.83        N/A        N/A        N/A
hydraulics            77.90        N/A        N/A        N/A
electric system      320.37        N/A        N/A        N/A
avionics             411.03        N/A        N/A        N/A
---------------------------------------------------------------
TOTAL               3625.03     4.1284     0.0000     0.0498
                                +39.89 % of MAC
---------------------------------------------------------------

Display all components

>>> ac.mass.display()
===============================================================
X45C total mass components breakdown
===============================================================
Name                |Mass,kg |CG:X,m    |CG:Y,m    |CG:Z,m    |
---------------------------------------------------------------
wing                 701.69     4.1333     0.0000     0.0000
fuselage             114.06     3.4550     0.0000     0.0000
main landing gear    277.61     5.0000     0.0000     0.0000
nose landing gear     48.99     1.0000     0.0000     0.0000
engine              1036.00     4.1150     0.0000     0.1000
engine mount          65.29     4.1150     0.0000     0.1000
engine section        15.67     4.1150     0.0000     0.1000
engine cooling        11.59     4.1150     0.0000     0.1000
air intake            95.69        N/A        N/A        N/A
oil cooling           17.15        N/A        N/A        N/A
engine controls        6.17     4.1150     0.0000     0.1000
engine starter        75.75        N/A        N/A        N/A
fuel system          253.23        N/A        N/A        N/A
flight controls       96.83        N/A        N/A        N/A
hydraulics            77.90        N/A        N/A        N/A
electric system      320.37        N/A        N/A        N/A
avionics             411.03        N/A        N/A        N/A
drop payload        1032.00     4.1150     0.0000     0.1200
Fuel                3000.00     4.0630     0.0000     0.0000
---------------------------------------------------------------
TOTAL               7657.03     4.0951     0.0000     0.0376
                                +39.10 % of MAC
---------------------------------------------------------------

For more details on mass and balance methods see related part of documentation


Aerodynamics, Stability and Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Methods for aerodynamic analysis: 

Drag: Aero09 (default) or Friction+Korn equation

Lift, induced drag: AVL

>>> velocity = 200.
>>> ac.get_drag(velocity)
array(0.007356930569642486)

Drag analysis method can be switched

>>> ac.set_drag_method(1) #set aero09
>>> Mach = 0.6
>>> ac.get_drag(Mach)
array(0.00737)

Switch method to Friction+Korn

>>> ac.set_drag_method(2) # set Friction+Korn
>>> ac.get_drag(Mach)
array(0.007468284411719359)

For lift and induced drag AVL is used. 

Trim analysis example is shown below. Note sometimes solution cannot be found 
due to improper CG, mass or elevator geometry input. Mass and CG may be 
specified by user or calculated by program.

>>> velocity = 200.
>>> altitude = 5000.
>>> results = ac.get_aero_trim(velocity, altitude)
>>> results.display()
AVL analysis results
====================
a        = +3.09092      
xNP      = +4.07722      
elevator = -5.08269      
e        = +0.76393      
ARref    = +2.92521      
Bref     = +9.74000      
k        = +0.14244      
Sref     = +32.43100     
beta     = +0.00000      
CD0      = +0.00741      
SM       = -0.00422      
alpha    = +0.74432      
CL0      = +0.03845      
Cref     = +4.24510      
Mach     = +0.62400      
Coefficients
------------
Cm       = +0.00000      
Cl       = +0.00000      
CL       = +0.07860      
CDind    = +0.00088      
CD       = +0.00830      
CZ       = -0.07870      
CY       = +0.00000      
CX       = -0.00727      
Cn       = +0.00000      
Derivatives
-----------
Clde     = +0.00000      
CYr      = -0.00046      
CYq      = +0.00000      
CYp      = +0.04129      
Cnr      = -0.00014      
Cnp      = -0.01081      
Cnq      = +0.00000      
Cla      = +0.00000      
Clb      = -0.03029      
CYb      = -0.00751      
CYa      = +0.00000      
Cnb      = +0.00037      
Cna      = +0.00000      
Clp      = -0.23989      
Clq      = +0.00000      
Clr      = +0.03181      
ede      = +0.17314      
CYde     = +0.00000      
Cmde     = -0.00707      
Cmb      = +0.00000      
Cma      = +0.01305      
Cnde     = +0.00000      
Cmr      = +0.00000      
Cmq      = -0.85060      
Cmp      = +0.00000      
CLde     = +0.02407      
CDfe     = +0.00034      
CDde     = +0.00000      
CLa      = +3.09092      
CLb      = +0.00000      
CLp      = +2.71952      
CLq      = +0.00000      
CLr      = +0.00000      


Perform aerodynamic analysis of given configuration (no trim)

>>> alpha = 2.0 # degrees
>>> beta = 1.0
>>> elevator = 3.0
>>> results = ac.get_aero_single_point(200,5000,alpha,beta,elevator)
>>> results.display()
AVL analysis results
====================
a        = +3.06572      
xNP      = +4.07139      
elevator = +3.00000      
e        = +1.04254      
ARref    = +2.92521      
Bref     = +9.74000      
k        = +0.10438      
Sref     = +32.43100     
beta     = +1.00000      
CD0      = +0.00741      
SM       = -0.00560      
alpha    = +2.00000      
CL0      = +0.23319      
Cref     = +4.24510      
Mach     = +0.62400      
Coefficients
------------
Cm       = -0.05679      
Cl       = -0.00161      
CL       = +0.34020      
CDind    = +0.01208      
CD       = +0.02005      
CZ       = -0.34069      
CY       = -0.00018      
CX       = -0.00816      
Cn       = +0.00002      
Derivatives
-----------
Clde     = -0.00011      
CYr      = -0.00235      
CYq      = -0.00818      
CYp      = +0.12202      
Cnr      = -0.00271      
Cnp      = -0.03961      
Cnq      = +0.00027      
Cla      = -0.01646      
Clb      = -0.08680      
CYb      = -0.00907      
CYa      = -0.00006      
Cnb      = +0.00377      
Cna      = +0.00211      
Clp      = -0.23834      
Clq      = +0.00000      
Clr      = +0.09069      
ede      = +0.00085      
CYde     = -0.00001      
Cmde     = -0.00705      
Cmb      = +0.00122      
Cma      = +0.01716      
Cnde     = +0.00000      
Cmr      = +0.00052      
Cmq      = -0.85092      
Cmp      = -0.00142      
CLde     = +0.02396      
CDfe     = +0.00169      
CDde     = +0.00000      
CLa      = +3.06572      
CLb      = -0.01001      
CLp      = +2.72113      
CLq      = -0.00013      
CLr      = -0.00097      


Propulsion
~~~~~~~~~~

wrapper function to jet engine analysis routine calculates TSFC values at 
given velocity (or Mach), altitude and thrust

>>> Mach = 0.5
>>> altitude = 1e4 # 10km
>>> requiredThrust = 2e4 # 20kN
>>> ac.get_sfc(Mach, altitude, requiredThrust)
0.00026576595734847475
"""
from db_tools import ReadDatabase
from paths import MyPaths
import numpy as np
from scipy.interpolate import Akima1DInterpolator

from engine_turbofan import Propulsion
from weight_fw import get_flying_wing_mass, AircraftMass
from flight_conditions import FlightConditions
from display_aircraft import flying_wing_display
from wing import Wing
from aero_avl_fw import Aerodynamics, FlightConditionsAVL
from drag_aero09 import get_parasite_drag_fw
from drag import get_transonic_drag_wing


import matplotlib.pyplot as plt

path = MyPaths()

def load(name):
    """
    Load aircraft from *.xls datasheet.
    
    Parameters
    ----------
    
    name : string
        aircraft name
    """
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
        self._dragAero09 = True # if True then aero09 is used, otherwise Mason+Korn equation

    def load_xls(self, name, xlsPath=None):
        """
        Loads aircraft data from *.xls datasheet and calculates necessary 
        parameters:
        - full geometry data set
        - weight and cg
        - parasite drag curve: CD0 vs Mach
        - engine thrust table if it does not exist for selected engine
        
        Parameters
        ----------
        
        name : string
            name of aircraft in database
        xlsPath : path, optional
            database path
        """
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
        self.propulsion.load(engineName,True)
        self.propulsion.sfcModel = None
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
        self._process_data(True)
#        fuelCG = self.wing.locate_on_wing(self.wing.fuelTankCGratio[0],self.wing.fuelTankCGratio[1])
#        fuelCG[1] = 0.0
#        self.mass.set_fuel_mass(self.designGoals.fuelMass,fuelCG)
        self.designGoals._process_data()
        self._update_parasite_drag()
        self._update_mass()

    def _process_data(self,updateAirfoils=False):
        self.wing._process_data(updateAirfoils)
        self.designGoals.set_flight_conditions(self.wing.MAC)
        self.designGoals._process_data()
        self.propulsion._process_data()

    def save_xls(self):
        """
        TBD.
        Saves configuration to *.xls datasheet
        """
        pass

    def display(self,showAxes=False):
        """
        Displays current aircraft configuration using mayavi. 
        
        Notes
        -----
        
        Method may produce an error depending on PC configuration. To fix the error 
        follow these steps:
        
        1. In tools->Preferences->Console->External modules
        set value of ETS_TOOLKIT to "wx"
        
        2. 	In C:/Python27/Lib/site-packages/tvtk/pyface/light_manager.py, 
        CameraLight class, _color_changed method, change: 
        
        >>> self.source.color = val
        
        to
        
        >>> self.color = val
        """
        flying_wing_display(self,showAxes)
    
    def display_2d(self,linetype='k-'):
        """
        Displays 2D planform shape of the wing.

        Parameters
        ----------
        
        linetype : string, optional
            matplotlib linetype format
        """
        plt.figure()
        plt.plot(self.wing.x2d, self.wing.y2d,linetype)
        plt.axis('equal')
        plt.show()
    
    def _update_parasite_drag(self):
        alt = self.designGoals.cruiseAltitude
        if self._dragAero09:
            M, CD, Mdd, CDdd = get_parasite_drag_fw(self,alt)
            self._M = np.hstack([M[0]-.2,M[0]-.1,M])+.1
            self._CD = np.hstack([CD[0],CD[0],CD])
            self.Mdd = Mdd
            self.CDdd = CDdd
        else:
            M, CD = get_transonic_drag_wing(self.wing,Ka=0.90)
            self._M = M
            self._CD = CD
        self._dragCurve = Akima1DInterpolator(self._M,self._CD)
    
    def plot_drag(self):
        m = np.linspace(self._M[0],self._M[-1],100)
        cd = np.array([self.get_drag(_m) for _m in m])
        plt.figure()
        plt.hold(True)
        #plt.plot(self._M, self._CD,'bs-')
        plt.plot(m,cd,'r-')
        plt.xlabel('Mach')
        plt.ylabel('Drag coefficient')
        plt.grid(True)
        plt.show()

    def _update_mass(self):
        self.mass = get_flying_wing_mass(self)
#    
#    def get_parasite_drag(self,velocity=None,altitude=None):
#        """
#        Returns friction and form drag coefficient at given velocity and altitude. 
#        If velocity and altitude are not specified then values from design targets
#        are used
#        
#        Parameters
#        ----------
#        
#        velocity : float, m/sec
#        altitude : float, m
#        """
#        self._update_parasite_drag()
#        return self.drag.get_total_drag()
    
    def update_aero_trim(self,velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                      mass=None,cg=None,inertia=None,CD0=None):
        """
        Updates results of aerodynamic analysis at trim condition using AVL solver. Stores the result 
        of analysis in self.aeroResults.
        
        Parameters
        ----------
        
        velocity : float, m/sec
            true airspeed of the aircraft. If velocity<5 then it is treated as Mach number, otherwise 
            airspeed in m/sec.
        altitude : float, m
            density altitude at which analysis is performed
        CmTrim : float
            Required moment coefficient. By default CmTrim=0 - no pitching moment.
        loadFactor : float
            load factor can be used for gust analysis: analysis is performed with mass= mass*loadFactor
        mass : float, kg
            aircraft mass. If value is not defined then total aircraft mass will be calculated.
        cg : array, m
            center of gravity in format array([x,y,z]). If value is not specified then value will be 
            calculated
        inertia : array, kg*m2
            aircraft moment of inertia in format array([Ixx, Iyy, Izz]). Required for dynamic 
            stability calculation.
        CD0 : float
            parasite drag coefficient. If value is not specified then value will be calculated.
        """
        aero = Aerodynamics(self)
        fc = FlightConditionsAVL(self,velocity,altitude,CmTrim,loadFactor,mass,
                                 cg,inertia,CD0)
        self.aeroResults = aero.run_trim(fc)
    
    def get_aero_trim(self,velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                      mass=None,cg=None,inertia=None,CD0=None):
        """
        Same as update_aero_trim, but returns results of analysis.
        """
        self.update_aero_trim(velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                              mass=None,cg=None,inertia=None,CD0=None)
        return self.aeroResults
    
    def get_aero_single_point(self,velocity,altitude,alpha=0.0,beta=0.0,
                              elevator=0.0,mass=None,cg=None,inertia=None,CD0=None):
        """
        Performs aerodynamic analysis using AVL at given aircraft configuration and flight condtions.
        
        
        Parameters
        ----------
        
        velocity : float, m/sec
            true airspeed of the aircraft. If velocity<5 then it is treated as Mach number, otherwise 
            airspeed in m/sec.
        altitude : float, m
            density altitude at which analysis is performed
        alpha : float, deg
            aircraft angle of attack
        beta : float, deg
            aircraft sideslip angle
        elevator : float, deg
            elevator deflection. positive direction is down.
        mass : float, kg
            aircraft mass. If value is not defined then total aircraft mass will be calculated.
        cg : array, m
            center of gravity in format array([x,y,z]). If value is not specified then value will be 
            calculated
        inertia : array, kg*m2
            aircraft moment of inertia in format array([Ixx, Iyy, Izz]). Required for dynamic 
            stability calculation.
        CD0 : float
            parasite drag coefficient. If value is not specified then value will be calculated.
        """
        aero = Aerodynamics(self)
        fc = FlightConditionsAVL(self,velocity,altitude,0,1,mass,cg,inertia,CD0)
        alpha = float(alpha)
        beta = float(beta)
        elevator = float(elevator)
        self.aeroResults = aero.run_single_point(fc,alpha,beta,elevator)
        return self.aeroResults

    def get_cg(self,update=True):
        """
        Returns center of gravity coordinates in format array([x,y,z]).
        
        Parameters
        ----------
        
        update : bool
            runs weight and balance analysis first. Use True if configuration has been changed since 
            last function call.
        """
        if update:
            self._update_mass()
        return self.mass.total.CG

    def get_mass(self,update=True):
        """
        Returns total mass (empty+payload) of current aircraft configuration.
        
        Parameters
        ----------
        
        update : bool
            runs weight and balance analysis first. Use True if configuration has been changed since 
            last function call.
        """
        if update:
            self._update_mass()
        return self.mass.total.totalMass

    def get_mass_empty(self,update=True):
        """
        Returns empty mass (airframe+propulsion) of current aircraft configuration.
        
        Parameters
        ----------
        
        update : bool
            runs weight and balance analysis first. Use True if configuration has been changed since 
            last function call.
        """
        if update:
            self._update_mass()
        return self.mass.empty.totalMass

    def set_drag_method(self,option=1):
        """
        Switches drag analysis methods.
        
        Parameters
        ----------
        
        option : int
            if 1 - Aero09, else Friction+Korn
        """
        if option==1:
            self._dragAero09 = True
        else:
            self._dragAero09 = False
        self._update_parasite_drag()

    def get_drag(self,velocity=None,altitude=None):
        """
        Returns parasite drag value of current aircraft configuration. Calculation is performed using 
        in-house Aero09 code.
        
        Parameters
        ----------
        
        velocity : float, m/sec
            true airspeed of the aircraft. If velocity<5 then it is treated as Mach number, otherwise 
            airspeed in m/sec.
        altitude : float, m
            density altitude at which analysis is performed
        """
        if velocity==None:
            velocity = self.designGoals.cruiseSpeed
        if altitude==None:
            altitude = self.designGoals.cruiseAltitude
        fc = FlightConditions(velocity,altitude)
        cd = self._dragCurve(fc.Mach)
        return cd
    
    def get_inertia(self):
        """
        Returns moment of inertia. Now this analysis is not available, so function returns [1;1;1].
        """
        return np.zeros(3) #FIXME: should be replaced by real calculation
    
    def get_sfc(self,velocity,altitude,thrustRequired):
        """
        Returns thrust specific fuel consumption (TSFC).
        
        Parameters
        ----------
        
        velocity : float, m/sec
            true airspeed of the aircraft. If velocity<5 then it is treated as Mach number, otherwise 
            airspeed in m/sec.
        altitude : float, m
            density altitude at which analysis is performed
        thrustRequired : float, N
            required thrust
        """
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


# --- debug and test section ---
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
    ac = load('Baseline1')
#    print ac.mass.empty()
#    ac.mass.display()
#    print ac.propulsion.get_sfc(0.7, 1e4, 2500)
    print ac.wing.get_max_segment_length(ac.propulsion.engine.diameter)
    print ac.propulsion.engine.length
    print ac.wing.segSweepTEdeg
    print ac.wing.segSweepLEdeg
    ac.mass.empty.display()
    ac.mass.display()
    #ac.display_2d()
    ac.plot_drag()
    ac.get_aero_single_point(0.7,1e4,0,0).display()
    print (ac.aeroResults.xNP - ac.wing.MAClocation)/ac.wing.MAC
    print ac.wing.MAClocation
    print ac.wing.MAC
    ac.display()

if __name__=="__main__":
    run_test4()