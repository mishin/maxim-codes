# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 11:38:33 2014

@author: Maxim
"""
from flight_conditions import FlightConditions, ISAtmosphere
import aircraft_FW

class BasicInput:
    """
    The class providing parameters required for most performance calculations. 
    Most parameters are calculated as properties of aircraft object.
    
    Parameters
    ----------
    
    aircraft : object
    Clmax : float
        maximum lift coefficient
    altitude : float, m
        altitude
    
    Arguments
    ---------
    
    density : float
        atmospheric density in kg/m**3
    Cd0 : float
        parasite drag constant
    k : float
        induced drag constant
    refArea : float
        wing reference area
    mass : float
        aircraft mass in kg
    g : float
        local gravity
    atm : ISAtmosphere
        ISA at given altitude
    Clmax : float
        maximum lift coefficient
    wt : float
        aircraft weight in N
    """
    def __init__(self,aircraft,Clmax=2.0,altitude=0.0):
        self.ac = aircraft
        self.tm = aircraft.propulsion
        self.atm = ISAtmosphere(altitude)
        self.g = self.atm.g
        self.density = self.atm.density
        self.Clmax = Clmax #FIXME: should be based on user's input or calculation
        self.mass = self.ac.mass()
        self.wt = self.ac.mass() * self.atm.g
        self.refArea = self.ac.wing.area
        self.update_aero()
    
    def set_CLmax(self,CLmaxNew):
        """
        sets new maximum lift coefficient
        """
        self.Clmax = float(CLmaxNew)

    def update_aero(self,velocity=None,altitude=None,flapDefl=0.0):
        """
        updates aerodynamic characteristics of given aircraft
        """
        if velocity==None:
            velocity = self.ac.designGoals.cruiseSpeed
        if altitude==None:
            altitude = self.ac.designGoals.cruiseAltitude
        self.atm = ISAtmosphere(altitude)
        self.density = self.atm.density
        self.Cd0 = self.ac.get_drag(velocity,altitude)
        self.ac.update_aero_trim(velocity,altitude)
        self.k = self.ac.aeroResults.k

    def update_fuel_mass(self,fuelMass):
        """
        Updates mass of fuel and CG due to this change
        """
        self.ac.mass.fuel.set_fuel_mass(fuelMass)
        self.mass = self.ac.mass()
        self.wt = self.mass*self.g
    
    def update_fuel_mass_burned(self,fuelMass):
        self.ac.mass.fuel.set_fuel_burned(fuelMass)
        self.mass = self.ac.mass()
        self.wt = self.mass*self.g

    def get_drag_coefficient(self,CL):
        """
        Calculates drag coefficient
        """
        return self.Cd0 + self.k*CL*CL

    def get_dynamic_pressure(self,V):
        """
        Calculates dynamic pressure
        """
        return self.atm.density*V*V/2.0




class FlightMechanics:
    def __init__(self,basicInput, thrustModule):
        self.bi = basicInput
        self.tm = thrustModule
    
    def get_Vstall(self,n=1.0):
        """
        Calculates stall speed at given load factor
        """
        Vs=(2.*self.bi.wt*n/self.bi.density/self.bi.refArea/self.bi.Clmax)**.5
        return Vs

    def get_EAS(self,V,rho):
        r"""
        Calculates equivalent airspeed
        .. math:
            V_{EAS} = V\sqrt{\frac{\rho}{\rho_{SL}}}        
        """
        atm0=ISAtmosphere(0.0)
        return V*(rho/atm0.density)**.5

    def get_required_thrust(self,V,altitude):
        r"""
        Determines required thrust at a given velocity and flight condition
        
        .. math:
            Lift = Weight
            q = \frac{\rho V^2}{2}
            C_L = \frac{Lift}{q S}
            C_D = C_{D0} + k C_L^2
            Drag = C_D S Q
            Thrust = Drag
        """
        atm = ISAtmosphere(altitude)
        rho = atm.density
        Cd0 = self.bi.ac.get_drag(V,altitude)
        g   = self.bi.g
        k   = self.bi.k
        S   = self.bi.refArea
        W   = self.bi.mass*g
        L   = W
        Q=0.5*rho*V*V
        Cl=L/Q/S
        Cd=Cd0+k*Cl**2
        D=Q*S*Cd
        return D

    def get_V_LDmax(self):
        r"""
        Calculates velocity at maximum lift-to-drag ratio
        """
        rho = self.bi.density
        W   = self.bi.wt
        S   = self.bi.refArea
        k   = self.bi.k
        Cd0=self.bi.Cd0
        return (2.*W/rho/S*(k/Cd0)**.5)**.5

    def get_LDmax(self):
        r"""
        Calculates maximum lift-to-drag ratio of an aircraft
        """
        return (1./4./self.bi.Cd0/ self.bi.k)**.5

# --- debug section ---

def run_test1():
    ac = aircraft_FW.load('X45C')
    bi = BasicInput(ac)
    fm = FlightMechanics(bi,ac.propulsion)
    print fm.get_Vstall()
    print fm.get_EAS(100.,1.0)
    print fm.get_LDmax()
    print fm.get_V_LDmax()
    print fm.get_required_thrust(100,1e4)

def run_test2():
    import matplotlib.pyplot as plt
    import numpy as np
    ac = aircraft_FW.load('X45C')
    bi = BasicInput(ac)
    fm = FlightMechanics(bi,ac.propulsion)
    V = np.linspace(100,400,50)
    Treq = np.array([fm.get_required_thrust(v,10000)/1e3 for v in V])
    plt.figure()
    plt.plot(V,Treq)
    plt.show()    

if __name__=="__main__":
    run_test2()