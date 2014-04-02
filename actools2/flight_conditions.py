# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:08:33 2014

@author: Maxim
"""
from math import exp
import constants as const

class FlightConditions:
    """
    stores flight conditions of an aircraft. Atmosphere properties are calculated
    using ISAtmosphere and stored in self.atmosphere
    
    Parameters
    ----------
    speed : float
        if speed > 3.0 then speed is treated as airspeed [m/sec]
        if speed < 3.0 then treated as Mach number
        
    altitude : float, meters
        altitude
    
    dT : degrees Kelvin
        temperature deviation from standard atmosphere

    refLength : float, meters
        reference length to calculate Reynolds number
    """
    def __init__(self,speed, altitude, dT=0.0, refLength=1.0):
        self.atm = ISAtmosphere(altitude,dT)
        self.refLength = float(refLength)
        if speed<=3.0:
            self.Mach = float(speed)
            self.velocity = self.Mach * self.atm.soundSpeed
        else:
            self.velocity = float(speed)
            self.Mach = self.velocity / self.atm.soundSpeed
        self.Re = self.atm.density*self.velocity*self.refLength/self.atm.viscosity
        self.dynamicPressure = self.atm.density*self.velocity*self.velocity/2.0
    
    def set_speed(self,speed):
        """
        sets new speed
        """
        self.__init__(float(speed),self.atm.altitude,self.atmosphere.dT,self.refLength)

    def set_altitude(self,altitude):
        """
        sets new altitude. Atmosphere parameters are recalculated using new altitude.
        """
        self.__init__(self.velocity,altitude,self.atm.dT,self.refLength)
        
    def set_ref_length(self,refLength):
        """
        sets new reference length
        """
        self.__init__(self.velocity,self.atm.altitude,self.atm.dT,float(refLength))
    
    def get_wall_spacing(self,yplus=1.0):
        """
        calculates minimum wall spacing based on yplus value.
        """
        Cf = 0.026 / (self.Re**(1/7))
        rho = self.atm.density
        tauWall = Cf*rho * self.velocity**2/2
        Ufric = (tauWall/rho)**0.5
        ds = yplus*self.atm.viscosity / (Ufric*rho)
        return ds


class ISAtmosphere:
    def __init__(self,altitude,dT,gas='air'):
        altitude = float(altitude)
        P0 = 101325.0
        T0 = 288.15
        G0 = 9.80665
        R  = 287.04
        if altitude<11000:
            T = T0 - 6.5 * altitude / 1000 + dT
            P = P0 * (1 - 0.0065 * altitude / T0) ** 5.2561
        else:
            T11 = T0 - 6.5 * 11000 / 1000 + dT
            P11 = P0 * (1 - 0.0065 * 11000 / T0) ** 5.2561
            T = T11 + dT
            P = P11 * exp(-G0 / (R * T11) * (altitude - 11000))  
        rho = P / (R * T)
        a = (1.4 * R * T) ** 0.5
        self.dT          = dT
        self.gas         = gas
        self.altitude    = altitude
        self.temperature = T
        self.pressure    = P
        self.soundSpeed  = a
        self.density     = rho
        self.set_sutherland_viscosity(self.gas)
        self.g = const.get_gravity_acceleration(altitude)
    
    def set_sutherland_viscosity(self,gas):
        if gas == 'ideal-gas':
            mu0 = 1.716e-5
            T0  = 273.11
            C   = 110.56
        elif gas == 'air':
            mu0 = 18.27e-6
            T0  = 291.15
            C   = 120.0
        mu = mu0*(T0+C)/(self.temperature+C)*((self.temperature/T0)**(3/2))
        self.viscosity = mu


def get_density_altitude(density):
    """
    Finds altitude with given density. 
    
    Parameters
    ----------
    density : float, kg/m^3
         atmospheric density
    
    Note
    ----
    approximated method is valid for altitude < 11km
    """
    a =  3.30858898e-09
    b = -1.13728614e-04
    c =  1.22202839e+00 - density
    return (-b - (b*b - 4.0*a*c)**0.5) / (2.0*a)


# --- debug ---
def run_test1():
    fc = FlightConditions(65.0,1500.)
    print fc.atm.density
    print fc.atm.viscosity
    print fc.atm.soundSpeed
    print fc.Re
    alt = get_density_altitude(1.0)
    fc.set_altitude(alt)
    print fc.atm.density
    
    print const.get_gravity_acceleration(10000)
if __name__=="__main__":
    run_test1()