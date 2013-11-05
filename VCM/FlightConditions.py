# -*- coding: utf-8 -*-
"""
Created on Thu Aug 09 21:45:22 2012
@author: Maxim
"""
import math

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
    def __init__(self, speed=15.0, altitude=0., dT=0, refLength=1.0):
        self.atmosphere = ISAtmosphere(altitude, dT)
        self.speed      = float(speed)
        if self.speed<=3.0:
            self.Mach     = float(speed)
            self.velocity = self.Mach * self.atmosphere.soundSpeed
        else:
            self.velocity = float(speed)
            self.Mach     = self.velocity / self.atmosphere.soundSpeed
        self.refLength = float(refLength)
        
        rho = self.atmosphere.density
        V   = self.velocity
        L   = self.refLength
        mu  = self.atmosphere.viscosity
        self.Re = rho * V * L / mu
    def set_speed(self,speed):
        """
        sets new speed
        """
        self.__init__(float(speed),self.atmosphere.altitude,self.atmosphere.dT,self.refLength)
    def set_altitude(self,altitude):
        """
        sets new altitude. Atmosphere parameters are recalculated using new altitude.
        """
        self.__init__(self.speed,altitude,self.atmosphere.dT,self.refLength)
    def set_ref_length(self,refLength):
        self.__init__(self.speed,self.atmosphere.altitude,self.atmosphere.dT,float(refLength))
    
    def get_wall_spacing(self,yplus=1.0):
        Cf = 0.026 / (self.Re**(1/7))
        rho = self.atmosphere.density
        tauWall = Cf*rho * self.velocity**2/2
        Ufric = (tauWall/rho)**0.5
        ds = yplus*self.atmosphere.viscosity / (Ufric*rho)
        return ds

class ISAtmosphere:
    def __init__(self, altitude, dT=0, gas='air'):
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
            P = P11 * math.exp(-G0 / (R * T11) * (altitude - 11000))  
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
        self.g=getGravityConstant(altitude)
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
    def set_altitude(self,altitude):
        self.__init__(altitude)
    def set_delta_ISA(self,dT):
        self.__init__(self.altitude,dT)

class ISAtmosphere_simple:
    def __init__(self,altitude,dT=0.0):
        P0 = 101325.0
        T0 = 288.15
        G0 = 9.80665
        R  = 287.04
        gamma = 1.4
        altitude = float(altitude)
        if altitude<11000:
            T = T0 - 6.5 * altitude / 1000.0 + dT
            P = P0 * (1 - 0.0065 * altitude / T0) ** 5.2561
        else:
            T11 = T0 - 6.5 * 11000 / 1000 + dT
            P11 = P0 * (1 - 0.0065 * 11000 / T0) ** 5.2561
            T = T11 + dT
            P = P11 * math.exp(-G0 / (R * T11) * (altitude - 11000))
        self.density     = P / (R*T)
        self.temperature = T
        self.pressure    = P
        self.soundSpeed  = (gamma*R*T)**0.5

def get_density_altitude(rho):
    """
    calculates approximate value of density altitude
    
    Returns
    -------
    altitude : m
    """
    a = 3.30858898e-09
    b = -1.13728614e-04
    c = 1.22202839e+00 - rho
    return (-b - (b*b - 4.0*a*c)**0.5) / (2.0*a)

def getDensityAltitude_exact(rho):
    """
    Solves the ISA equations to determine the equivalent density
    
    Parameters
    ----------
    rho : float, kg/m^3
         atmospheric density 
    """
    def _f(h):
        atm=ISAtmosphere_simple(h)
        return atm.density-rho
    from scipy import optimize
    alt=optimize.broyden1(_f,0.0)
    return alt

def getGravityConstant(altitude):
    R=6371.0*10**3
    g=9.807
    c=R**2/(R+altitude)**2
    return g*c

# --- test functions ---
def run_test():
    alt = [0.0,2000.0,5000.0,7500.0,10000.0,11000.0]
    for h in alt:
        fc = FlightConditions(50.0,h)
        g=getGravityConstant(h)
        print 'Re = %.2e\tM = %.2f'%(fc.Re,fc.Mach)
        print 'g  = %.4f'%(g)

def run_test2():
    from numpy import linspace, polyfit
    import matplotlib.pyplot as plt
    alt = linspace(0,11000.0)
    density = list()
    for h in alt:
        fc = ISAtmosphere_simple(h)
        density.append(fc.density)
    n = 2
    coef = polyfit(alt,density,n)
    densityFit = list()
    for h in alt:
        _density = 0.0
        for i in range(n+1):
            _density += h**(n-i) * coef[i]
        densityFit.append(_density)
    print coef
    
    rho = 0.5
    alt1 = get_density_altitude(rho)
    
    plt.figure(1)
    plt.plot(alt,density,'*')
    plt.hold(True)
    plt.plot(alt, densityFit,'r--')
    plt.plot(alt1,rho,'go')
    plt.show()

if __name__=="__main__":
    run_test2()