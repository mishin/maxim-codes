# -*- coding: utf-8 -*-
"""
Contains classes for aircraft performance analysis
"""
import FlightConditions as fc
from scipy.optimize import minimize_scalar
import numpy as np
import aircraft
from miscTools import Timer

timer = Timer()

class performance:
    """
    incomplete
    
    This class wraps aircraft performance calculations for easy access
    Sets and runs default aerodyamics calculations
    
    Arguments
    ---------
    
    aircraft : aircraft object
    
    aero : aerodynamics.results
        results of aerodynamic analysis
    basicInput : basic input object
        basic input object at test flight conditions
    result : result
        result of requested analysis
    """
    def __init__(self, aircraft):
        self._atm=fc.ISAtmosphere(aircraft.designGoals.designAltitude)
        self.aircraft=aircraft
        self.aero=None
        self.basicInput=None
        self.result=None
    def setDefaultAerodynamics(self):
        self.setAerodynamics(self.aircraft.designGoals.designSpeed,self.aircraft.designGoals.designAltitude,0.0)
    def setAerodynamics(self,velocity,altitude,deltaISA,mass=None,flapSetting=0.0):
        if mass==None:
            mass=self.aircraft.get_mass_total()
        atm=fc.ISAtmosphere(altitude,deltaISA)
        self.aircraft.analysis.aerodynamics.update(velocity,atm.density)
        self.aero=self.aircraft.analysis.aerodynamics.results
        Cd0=self.aero.Cd0
        k=self.aero.k
        S=self.aircraft.wing.area
        Clmax=self.aircraft.analysis.aerodynamics.Clmax(velocity/atm.soundSpeed,flapSetting)
        self.basicInput=basicInput(altitude,deltaISA,Cd0,k,S,mass,Clmax)
    def getVelocity(self,powerSetting):
        if self.aero==None or self.basicInput==None:self.setDefaultAerodynamics()
        self.result=steadyLevelFlight()
        self.result.run_velocity_maxTAS(self.basicInput,self.aircraft.analysis.thrust,powerSetting)

class _basicInput:
    """
    Module will be removed
    This class holds common parameters required for running most performance 
    calculations

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
    
    Parameters
    ----------
    
    altitude : float
        altitude in m
    dT : float
        temperature deviation from standard ISA in degree C(K)
    Cd0 : float
        parasite drag constant
    k : float
        induced drag constant
    S : float
        wing reference area in m**2
    mass : float
        aircraft mass
    Clmax : float
        maximum lift coefficient
    """
    def __init__(self,altitude,dT,Cd0,k,S,mass,Clmax):
        atm=fc.ISAtmosphere(altitude,dT)
        self.density =atm.density
        self.Cd0     =Cd0
        self.k       =k
        self.refArea =S
        self.mass    =mass
        self.g       =atm.g
        self.atm     =atm
        self.Clmax   =Clmax
        self.wt      =self.mass*self.g
    def get_drag_coefficient(self,CL):
        """
        Returns total drag coefficient at given lift as CD = CD0 + k*CL^2
        """
        return self.Cd0 + self.k*CL*CL
        
        
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
        self.atm = fc.ISAtmosphere(altitude)
        self.g = self.atm.g
        self.density = self.atm.density
        self.Clmax = Clmax
        self.mass = self.ac.mass()
        self.wt = self.ac.mass() * self.atm.g
        self.refArea = self.ac.wing.area
        self.update_aero()
    
    def set_CLmax(self,CLmaxNew):
        """
        sets new maximum lift coefficient
        """
        self.Clmax = float(CLmaxNew)

    def update_aero(self,velocity=0.0,altitude=0.0,flapDefl=0.0):
        """
        updates aerodynamic characteristics of given aircraft
        """
        self.atm = fc.ISAtmosphere(altitude)
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
    """
    This class contains common functions for flight mechanics calculations
    and attributes for algorithm output
    
    Arguments
    ---------
    
    velocity : float
        true airspeed in m/sec
    equivalentAirspeed : float
        equivalent airspeed in m/sec
    densityAltitude : float
        density altitude in m
    density : float
        air density in kg/m**3
    climbRate : float
        climb rate in m/sec
    climbAngle : float
        climb angle in degrees
    fuelFlow : float
        fuel flow in kg/sec
    RPM : float
        engine revolutions per minute
    powerSetting : float [0;100]
        power setting in percent power
    thrust : float
        thrust in N
    power : float
        power in Watt
    SAR_km_per_kg : float
        specific air range in km/kg
    lift : float
        lift in N
    drag : float
        drag in N
    L_D : float
        lift to drag ratio
    turnRadius : float
        turn radius in meters
    bankAngle : float
        bank angle in radians
    loadFactor : float
        load factor in number of multiplies of earth's gravity
    turnRate : float
        turn rate in rad/sec
    """
    def __init__(self, basicInput, thrustModule):
        self.bi = basicInput
        self.tm = thrustModule
        self.velocity          =0.0
        self.equivalentAirspeed=0.0
        self.densityAltitude   =0.0
        self.density           =0.0
        self.climbRate         =0.0
        self.climbAngle        =0.0
        self.fuelFlow          =0.0
        self.RPM               =0.0   
        self.powerSetting      =0.0   
        self.thrust            =0.0 
        self.power             =0.0 
        self.propEfficiency    =0.0 
        self.beta              =0.0 
        self.SAR_km_per_kg     =0.0
        self.lift              =0.0
        self.drag              =0.0
        self.L_D               =0.0   
        self.turnRadius        =0.0
        self.bankAngle         =0.0
        self.loadFactor        =1.0
        self.turnRate          =0.0
    def vStall(self,n=1.0):
        """
        Calculates stall speed at given load factor
        """
        Vs=(2.*self.bi.wt*n/self.bi.density/self.bi.refArea/self.bi.Clmax)**.5
        return Vs
    def get_Vstall(self,n=1.0):
        return self.vStall(n)

    def EAS(self,V,rho):
        r"""
        Calculates equivalent airspeed
        .. math:
            V_{EAS} = V\sqrt{\frac{\rho}{\rho_{SL}}}        
        """
        atm0=fc.ISAtmosphere(0.0)        
        return V*(rho/atm0.density)**.5
    def get_EAS(self,V,density):
        return self.EAS(V,density)

    def requiredThrust(self,V):
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
        rho = self.bi.density
        Cd0 = self.bi.Cd0
        g   = self.bi.g
        k   = self.bi.k
        S   = self.bi.refArea
        W   = self.bi.mass*g
        L   = W
        Q=0.5*rho*V**2
        Cl=L/Q/S
        Cd=Cd0+k*Cl**2
        D=Q*S*Cd
        return D
    def get_requiredThrust(self,V):
        return self.requiredThrust(V)

    def V_LDmax(self):
        r"""
        Calculates velocity at maximum lift-to-drag ratio
        """
        rho = self.bi.density
        W   = self.bi.wt
        S   = self.bi.refArea
        k   = self.bi.k
        Cd0=self.bi.Cd0
        return (2.*W/rho/S*(k/Cd0)**.5)**.5
    def get_V_LDmax(self):
        return self.V_LDmax()

    def LDmax(self):
        r"""
        Calculates maximum lift-to-drag ratio of an aircraft
        """
        return (1./4./self.bi.Cd0/ self.bi.k)**.5
    def get_LDmax(self):
        return self.LDmax()

    def display(self):
        for attr, value in self.__dict__.iteritems():
            print '{:<20} = {:<12}  '.format(attr,str(value))


class ClimbAndDescendingFlight(FlightMechanics):
    def run_maximumClimbRate(self,powerSetting):
        """
        Calculates the maximum climb rate at a given throttle setting by 
        varying the airspeed
        """
        vStall                 =self.vStall()
        bound                  =np.array([vStall,2.0*vStall])
        Args                   =(powerSetting,)
        opts                   ={'maxiter':100,'disp':False}
        fun                    =self._clmEqn_maxRate
        res                    =minimize_scalar(fun, tol=0.01,bracket=bound, bounds=bound, args=Args, method='Bounded',options=opts)
        airspeed               =res.x      
        self.run_climbRate(airspeed,powerSetting)
    def run_maximumClimbAngle(self,powerSetting):
        """
        Calculates the maximum climb angle at a given throttle setting by 
        varying the airspeed
        """
        vStall                 =self.vStall()
        bound                  =np.array([vStall,2.0*vStall])
        Args                   =(powerSetting,)
        opts                   ={'maxiter':100,'disp':False}
        fun                    =self._clmEqn_maxAngle
        res                    =minimize_scalar(fun, bracket=bound, bounds=bound, args=Args, method='Bounded',options=opts)
        airspeed               =res.x      
        self.run_climbRate(airspeed,powerSetting)
    def run_minGlideAngle(self):
        """
        Calculates the minimum glide angle by varying the airspeed.  Useful
        for maximum glide range calculations
        """
        vStall                 =self.vStall()
        bound                  =np.array([vStall,2.0*vStall])
        opts                   ={'maxiter':100,'disp':False}
        fun                    =self._glide_minAngle
        res                    =minimize_scalar(fun, bracket=bound, bounds=bound, method='Bounded',options=opts)
        airspeed               =res.x
        self.run_climbRate(airspeed,0.)
    def run_minGlideSinkRate(self):
        """
        Calculates the minimum glide angle by varying the airspeed.  Useful
        for maximum glide endurance calculations
        """
        vStall                 =self.vStall()
        bound                  =np.array([vStall,2.0*vStall])
        opts                   ={'maxiter':100,'disp':False}
        fun                    =self._glide_minSinkRate
        res                    =minimize_scalar(fun, bracket=bound, bounds=bound, method='Bounded',options=opts)
        airspeed               =res.x      
        self.run_climbRate(airspeed,0.)
    def run_climbRate(self,airspeed,powerSetting):
        """
        Calculates the climb rate at a given airspeed and power setting.
        Supports descent calculations which will appear as negative climb
        rates. Does not utilize small angle approximations.
        """
        basicInput = self.bi
        thrustModule = self.tm
        V       =airspeed
        rho     =basicInput.density
        Cd0     =basicInput.Cd0
        g       =basicInput.g
        k       =basicInput.k
        S       =basicInput.refArea
        W       =basicInput.mass*g
        Q       =0.5*rho*V*V
        tol     =1.0e-4
        change  =tol+1.0
        CA      =0.0 #climb angle, rad
        iMax    =100
        i       =0
        if powerSetting>=5.0:
            tm=thrustModule.runAnalysis
        else:
            tm=thrustModule._zeroThrustModule
        while change>tol and i<=iMax:            
            tm(V,rho,powerSetting)
            T  =thrustModule.thrust
            L  =W*np.cos(CA)-T*np.sin(CA)
            Cl =L/Q/S
            Cd =Cd0+k*Cl*Cl
            D  =Q*S*Cd
            PR =D*V
            PA =T*V
            RC =(PA-PR)/W
            ca =CA
            CA =np.arcsin(RC/V)
            change=np.abs(ca-CA)
            i+=1
        self.__init__(self.bi, self.tm)
        self.velocity          =V
        self.equivalentAirspeed=self.EAS(V,basicInput.density)
        self.densityAltitude   =fc.get_density_altitude(basicInput.density)
        self.density           =basicInput.density
        self.climbRate         =RC
        self.climbAngle        =CA
        self.fuelFlow          =thrustModule.fuelFlow
        self.RPM               =thrustModule.RPM
        self.powerSetting      =powerSetting
        self.thrust            =T
        self.power             =PA
        self.propEfficiency    =thrustModule.propEfficiency
        self.beta              =thrustModule.beta
        self.SAR_km_per_kg     =thrustModule.SAR_km_per_kg
        self.lift              =L
        self.drag              =D
        self.L_D               =L/D   
        self.turnRadius        =0.0
        self.bankAngle         =0.0
        self.loadFactor        =1.0
    def _clmEqn_maxRate(self,airspeed,powerSetting):
        self.run_climbRate(airspeed,powerSetting)
        return -self.climbRate
    def _clmEqn_maxAngle(self,airspeed,powerSetting):
        self.run_climbRate(airspeed,powerSetting)
        return -100.*self.climbAngle
    def _glide_minAngle(self,airspeed):
        self.run_climbRate(airspeed,0.)
        return np.abs(self.climbAngle) 
    def _glide_minSinkRate(self,airspeed):
        self.run_climbRate(airspeed,0.)
        return np.abs(self.climbRate)


class SteadyLevelFlight(FlightMechanics):
    def run_velocity_maxTAS(self,powerSetting,bound=[20.,100.]):
        """
        Calculates the maximum aircraft velocity attainable at a given 
        flight condition and throttle setting.  Finds the propeller
        configuration and engine RPM that maximizes velocity
        """
        Args                   =(powerSetting,)
        opts                   ={'maxiter':100,'disp':False}
        fun                    =self._velocityEqn_maxTAS
        res                    =minimize_scalar(fun, bracket=bound, bounds=bound, args=Args, method='Bounded',options=opts)
        velocity               =res.x
        self._setAttributes(velocity)

    def run_velocity_maxSAR(self):
        """
        Calculates the velocity and throttle setting that maximizes the 
        distance covered per unit of fuel burned (km/kg).  This value can be 
        used as a part of maximum range calculations. The result should be a 
        better estimate of the true best range speed than (L/D)max since 
        propeller efficiency is taken into account
        """
        vStall                 =self.vStall()
        bound                  =np.array([vStall,1.5*vStall])
        opts                   ={'maxiter':100,'disp':False}
        fun                    =self._velocityEqn_maxSAR
        res                    =minimize_scalar(fun, bracket=bound, bounds=bound, method='Bounded',options=opts)
        velocity               =res.x
        self._setAttributes(velocity)

    def run_velocity_minFuel(self):
        """
        Calculates the velocity that minimizes the fuel flow for level
        flight.  This value can be used as a part of maximum endurance 
        calculations.  The result should be a better estimate of the true
        best endurance speed than (CL^3/2 / CD)max since propeller efficiency
        is taken into account
        """
        vStall                 =self.vStall()
        bound                  =np.array([vStall,1.5*vStall])
        opts                   ={'maxiter':100,'disp':False}
        fun                    =self._velocityEqn_minFuel
        res                    =minimize_scalar(fun, bracket=bound, bounds=bound, method='Bounded',options=opts)
        velocity               =res.x      
        self._setAttributes(velocity)

    def _velocityEqn_maxTAS(self,V,powerSetting):
        self.tm.runAnalysis(V,self.bi.density,powerSetting)
        T=self.tm.thrust
        D=self.requiredThrust(V)
        return (T-D)**2
    def _velocityEqn_maxSAR(self,V):
        D=self.requiredThrust(V)
        self.tm.matchRequiredThrust(V,self.bi.density,D)
        FF=self.tm.fuelFlow
        SAR=V/FF
        return -SAR
    def _velocityEqn_minFuel(self,V):
        D=self.requiredThrust(V)
        self.tm.matchRequiredThrust(V,self.bi.density,D)
        FF=self.tm.fuelFlow
        return 10000.*FF
    def _setAttributes(self,velocity):
        drag=self.requiredThrust(velocity)
        self.__init__(self.bi, self.tm)
        self.tm.matchRequiredThrust(velocity,self.bi.density,drag)
        self.velocity          =velocity
        self.equivalentAirspeed=self.get_EAS(velocity,self.bi.density)
        self.densityAltitude   =fc.get_density_altitude(self.bi.density)
        self.density           =self.bi.density
        self.fuelFlow          =self.tm.fuelFlow
        self.RPM               =self.tm.RPM
        self.powerSetting      =self.tm.powerSetting
        self.thrust            =self.tm.thrust
        self.power             =self.tm.power
        self.propEfficiency    =self.tm.propEfficiency
        self.beta              =self.tm.beta
        self.SAR_km_per_kg     =self.tm.SAR_km_per_kg
        self.drag              =drag
        self.lift              =self.bi.wt
        self.L_D               =self.lift/self.drag


class TurningFlight(FlightMechanics):
    def run_turnRate(self,phi,airspeed,nMax,powerSetting):
        """
        Calculates the turn rate at a given airspeed, bank angle, and power 
        setting.
        """     
        nReq                   =1./np.cos(phi)
        nMax                   =self._limitLoadFactor(airspeed,nMax,powerSetting)        
        n                      =np.min([nMax,nReq])        
        self._setAttributes(n,airspeed,powerSetting)
                      
    def _limitLoadFactor(self,airspeed,nMax,powerSetting):
        V         =airspeed
        rho       =self.bi.density
        Cd0       =self.bi.Cd0
        g         =self.bi.g
        k         =self.bi.k
        S         =self.bi.refArea
        W         =self.bi.mass*g        
        Clmax     =self.bi.Clmax
        Q         =0.5*rho*V**2
        
        self.tm.runAnalysis(V,rho,powerSetting)
        T         =self.tm.thrust                
        # Thrust Limited Turn
        nThrust   =(Q/(k*W/S)*T/W-Q*Cd0/(W/S))**.5
        # Stall Limited Turn
        nStall    =Q*Clmax/(W/S)
        # Maximum Turn Rate Load Factor
        return np.min([nThrust,nStall,nMax])    
    def _setAttributes(self,n,airspeed,powerSetting):
        V                      =airspeed
        rho                    =self.bi.density
        Cd0                    =self.bi.Cd0
        g                      =self.bi.g
        k                      =self.bi.k
        S                      =self.bi.refArea
        W                      =self.bi.mass*g        
        Q                      =0.5*rho*V**2
        self.tm.runAnalysis(V,rho,powerSetting)
        phi                    =np.arccos(1./n)
        turnRadius             =V**2/(g*(n**2-1.)**.5)
        turnRate               =g*(n**2-1)**.5/V
        L                      =n*W
        Cl                     =L/Q/S
        Cd                     =Cd0+k*Cl**2
        D                      =Q*S*Cd
        self.tm.matchRequiredThrust(V,self.bi.density,D)
        self.__init__(self.bi,self.tm)
        self.velocity          =V
        self.equivalentAirspeed=self.EAS(V,self.bi.density)
        self.densityAltitude   =fc.get_density_altitude(self.bi.density)
        self.density           =self.bi.density
        self.fuelFlow          =self.tm.fuelFlow
        self.RPM               =self.tm.RPM   
        self.powerSetting      =powerSetting 
        self.thrust            =self.tm.thrust
        self.power             =self.tm.power
        self.propEfficiency    =self.tm.propEfficiency
        self.beta              =self.tm.beta
        self.SAR_km_per_kg     =self.tm.SAR_km_per_kg
        self.lift              =L
        self.drag              =D
        self.L_D               =L/D   
        self.turnRadius        =turnRadius
        self.bankAngle         =phi
        self.loadFactor        =n
        self.turnRate          =turnRate


def takeoff_rotation(ac,fieldAlt=0.0,stallFrac=0.8,mu=0.02,flapDefl=0.0):
    """
    Calculates required elevator deflection for takeoff rotation. Moment required
    for takeoff rotation is calculated as sum of the mass and friction moment 
    around main landing gear contact point. Returns elevator deflection angle 
    required to rotate aircraft at takeoff.
    
    Parameters
    ----------

    ac : aircraft
        aircraft object
    fieldAlt : float, m
        field (taxiway) altitude in meters
    stallFrac : float [0,1]
        fraction of stall velocity at which takeoff rotation is estimated
    mu : float
        ground friction coefficient
    """
    #FIXME: assumed values, replace with calculations, while loop can be added
    Clmax = 2.2
    ClTO = 0.9*Clmax
    CG = [1.5662, 0, 0]

    S = ac.wing.area
    print S
    bIn = BasicInput(ac,Clmax,fieldAlt)
    fm = FlightMechanics(bIn,ac.analysis.thrust)
    Vstall = fm.vStall()
    Vstall = 23.15
    V = Vstall*stallFrac
    q = bIn.get_dynamic_pressure(V)
    MLGx = ac.landingGear.groundContact_X[1]
    MLGz = ac.landingGear.groundContact_Z[1]
    
    if CG[0]>MLGx:
        raise ValueError('MLG should be located after CG')
    lift = 0.5*bIn.atm.density*V*V*bIn.refArea*ClTO
    momentMass = (MLGx-CG[0])*bIn.wt
    momentFrict = (CG[2]-MLGz)*mu*(bIn.wt-lift)
    CmRequired = (momentMass+momentFrict) / (q*S*ac.wing.MAC)
    print CmRequired, bIn.wt
    aeroRslt = ac.analyze_aero_trim(V,fieldAlt,1.0,-CmRequired,flapDefl)
    elevDeflection = aeroRslt.elevator
    return elevDeflection


def run_takeoff_rotation():
    ac = aircraft.load('V200')
    print takeoff_rotation(ac)

def run_performance_test(inputSheet='V0510'):
    from report_tools import header
    velocity =50.
    altitude =2000.
    Clmax    =1.4
    nMax     =4.5
    
    ac       =aircraft.load(inputSheet)

    tm       =ac.analysis.thrust
    bi       =BasicInput(ac, Clmax, altitude)
    
    slf      = SteadyLevelFlight(bi,tm)
    clm      = ClimbAndDescendingFlight(bi,tm)
    trn      = TurningFlight(bi,tm)
    
    print header('TURN')
    trn.run_turnRate(60.*np.pi/180.,velocity,nMax,100)
    trn.display()
    print header('CLIMB RATE')
    clm.run_climbRate(velocity,100)
    clm.display()
    print header('MAXIMUM CLIMB ANGLE')
    clm.run_maximumClimbAngle(100.)
    clm.display()
    timer.start()
    print header('MAXIMUM CLIMB RATE')
    clm.run_maximumClimbRate(100.)
    clm.display()
    timer.finish()
    print header('MINIMUM GLIDE ANGLE')
    clm.run_minGlideAngle()
    clm.display()
    print header('MINIMUM GLIDE SINK RATE')
    clm.run_minGlideSinkRate()
    clm.display()
    
    print header('MAXIMUM VELOCITY')
    slf.run_velocity_maxTAS(75.)
    slf.display()
    timer.start()
    print header('MAXIMUM SAR')
    slf.run_velocity_maxSAR()
    slf.display()
    timer.finish()
    print header('MINIMUM FUEL')
    slf.run_velocity_minFuel()
    slf.display()

def performanceTest():
    ac       =aircraft.load('sampleInput1')
    print 'Performance Wrapper Test'
    ac.analysis.performance.getVelocity(100.)
    ac.analysis.performance.result.display()
    print '='*20

def debug1():
    ac = aircraft.load('V0510')
    print takeoff_rotation(ac,0.0)

def debug2():
    ac = aircraft.load('V0510')
    extInp = BasicInput(ac,1500.0,2.1)
    print extInp

def test_flight_mechanics():
    ac = aircraft.load('V0510')
    bi = BasicInput(ac,1.6,2000.)
    tm = ac.analysis.thrust
    fm = FlightMechanics(bi,tm)
    
    print fm.get_EAS(50,1.05)
    print fm.get_LDmax()
    print fm.get_requiredThrust(65)
    print fm.get_Vstall(1.2)

if __name__=='__main__':
    #test_flight_mechanics()
    run_performance_test()
    