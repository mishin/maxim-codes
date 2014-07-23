# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 11:44:42 2014

@author: Maxim
"""

from perf_tools import FlightMechanics, ISAtmosphere
from scipy.optimize import fsolve, minimize_scalar
import numpy as np
import constants

class PerformanceResults:
    def __init__(self):
        self.velocity        = 0.0
        self.Mach            = 0.0
        self.SAR             = 0.0
        self.EAS             = 0.0
        self.altitude        = 0.0
        self.density         = 0.0
        self.climbRate       = 0.0
        self.climbAngle      = 0.0
        self.fuelFlow        = 0.0
        self.TSFC            = 0.0
        #self.powerSetting    = 0.0
        self.thrust          = 0.0
        #self.power           = 0.0
        self.lift            = 0.0
        self.drag            = 0.0
        self.LD              = 0.0
        self.turnRadius      = 0.0
        self.bankAngle       = 0.0
        self.loadFactor      = 1.0
        self.turnRate        = 0.0

    def __repr__(self):
        #out = 'Aircraft performance\n====================\n'
        out = 'performance data\n=========\n'
        for attr, value in self.__dict__.iteritems():
            if type(value) is float or type(value) is np.float64:
                out += '{:<15} = {:<12.5f}\n'.format(attr,value)
        return out

    def display(self):
        print self.__repr__()



class SteadyLevelFlight(FlightMechanics):
    def run_max_TAS(self,altitude):
        """ Maximum velocity """
        atm = ISAtmosphere(altitude)
        def velocity_eqn(V,altitude):
            Treq = self.get_required_thrust(V,altitude)
            Tav = self.tm.get_thrust_available(altitude)
            return Treq-Tav
        V0 = atm.soundSpeed*0.8
        V = fsolve(velocity_eqn,V0,(altitude,))[0]
        self.set_results(V,altitude)
        self.results.thrust = self.tm.totalThrust
        return self.results
    
    def run_min_TAS(self,altitude):
        atm = ISAtmosphere(altitude)
        def velocity_eqn(V,altitude):
            Treq = self.get_required_thrust(V,altitude)
            Tav = self.tm.get_thrust_available(altitude)
            return Treq-Tav
        V0 = atm.soundSpeed*0.3
        V = fsolve(velocity_eqn,V0,(altitude,))[0]
        self.set_results(V,altitude)
        self.results.thrust = self.tm.totalThrust
        return self.results
    
    def run_max_SAR(self,altitude):
        """ Maximum range """
        atm = ISAtmosphere(altitude)
        Vstall = self.get_Vstall()
        args = (altitude,atm.soundSpeed,)
        bound = np.array([Vstall,0.85*atm.soundSpeed])
        opts = {'maxiter':100,'disp':False}
        def velocity_eqn(V,altitude,soundSpeed):
            D = self.get_required_thrust(V,altitude)
            sfc = self.tm.get_sfc(V/soundSpeed,altitude,D)
            ff = D*sfc #fuel flow kg/sec
            return -V/ff
        rslt = minimize_scalar(velocity_eqn,bracket=bound,bounds=bound,method='Bounded',options=opts,args=args)
        V = rslt.x
        self.set_results(V,altitude)
        return self.results
    
    def run_min_fuel(self,altitude):
        """ Maximum endurance """
        atm = ISAtmosphere(altitude)
        Vstall = self.get_Vstall()
        args = (altitude,atm.soundSpeed,)
        bound = np.array([Vstall,0.85*atm.soundSpeed])
        opts = {'maxiter':100,'disp':False}
        def velocity_eqn(V,altitude,soundSpeed):
            D = self.get_required_thrust(V,altitude)
            sfc = self.tm.get_sfc(V/soundSpeed,altitude,D)
            ff = sfc*D
            return ff
        rslt = minimize_scalar(velocity_eqn,bracket=bound,bounds=bound,method='Bounded',options=opts,args=args)
        V = rslt.x
        self.set_results(V,altitude) #FIXME: results are not satisfactory
        return self.results
    
    def set_results(self,velocity,altitude):
        results = PerformanceResults()
        atm = ISAtmosphere(altitude)
        results.velocity = velocity
        results.Mach = velocity/atm.soundSpeed
        results.altitude = altitude
        results.density = atm.density
        results.thrust = self.get_required_thrust(velocity,altitude)
        results.drag = results.thrust
        results.TSFC = self.tm.get_sfc(results.Mach,altitude,results.thrust)
        results.fuelFlow = results.TSFC*results.thrust/constants.GRAVITY_ACCEL
        results.SAR = velocity/results.fuelFlow
        results.lift = self.bi.wt
        results.LD = results.lift / results.drag
        self.results = results


class ClimbDescent(FlightMechanics):
    def run_max_climb_rate(self,altitude):
        atm = ISAtmosphere(altitude)
        """ calculates maximum climb rate """
        Vstall = self.get_Vstall()
        bound = np.array([Vstall,0.9*atm.soundSpeed])
        opts = {'maxiter':100,'disp':False}
        args = (altitude,)
        fun = lambda velocity,altitude: -self.run_climb_rate(velocity,altitude).climbRate
        rslt = minimize_scalar(fun,bracket=bound,bounds=bound,method='Bounded',
                               options=opts,args=args)
        velocity = rslt.x
        return self.run_climb_rate(velocity,altitude)

    
    def get_service_ceiling(self):
        RCmin = 0.5
        fun = lambda h: self.run_max_climb_rate(h).climbRate-RCmin
        xopt = fsolve(fun,5e3,xtol=0.1)[0]
        return xopt

    def run_max_climb_angle(self):
        pass

    def run_min_glide_angle(self):
        pass
    def run_min_glide_sinkrate(self):
        pass
    def run_climb_rate(self,airspeed,altitude):
        basicInput   = self.bi
        thrustModule = self.tm
        atm = ISAtmosphere(altitude)
        V       =airspeed
        rho     =atm.density
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
        while change>tol and i<=iMax:            
            T = thrustModule.get_thrust_available(altitude)
            L  =W*np.cos(CA)-T*np.sin(CA)
            Cl =L/Q/S
            Cd =Cd0+k*Cl*Cl
            D  =Q*S*Cd
            RC =(T-D)*V/W
            ca =CA
            CA =np.arcsin(RC/V)
            change=np.abs(ca-CA)
            i+=1
        results = PerformanceResults()
        results.altitude = altitude
        results.velocity = airspeed
        results.thrust = T
        results.climbAngle = CA
        results.climbRate = RC
        results.Mach = airspeed/atm.soundSpeed
        results.density = atm.density
        results.lift = L
        results.drag = D
        results.LD = L/D
        results.TSFC = self.tm.get_sfc(results.Mach,altitude,results.thrust)
        results.fuelFlow = results.TSFC*results.thrust
        return results


class TurningFlight:
    def run_turn_rate(self,phi,airspeed,loadFactorMax,powerSetting):
        pass
    

class Field:
    pass

def run_test1():
    import aircraft_FW
    ac = aircraft_FW.load('Baseline1')
    ac.display()
    slf =SteadyLevelFlight(ac)
    #alt = ac.designGoals.cruiseAltitude
    alt = 10000.

    print 'MAXIMUM SPEED'
    print slf.run_max_TAS(alt)
    print 'MINIMUM SPEED'
    print slf.run_min_TAS(alt)
    print 'MAXIMUM SAR'
    print slf.run_max_SAR(alt)
    print 'MINIMUM FUEL'
    print slf.run_min_fuel(alt)

    clm = ClimbDescent(ac)
    #print clm.run_max_climb_rate(alt)
    print clm.get_service_ceiling()
    

if __name__=="__main__":
    run_test1()