# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 16:50:41 2014

@author: Maxim
"""

from performance import SteadyLevelFlight, ClimbDescent, BasicInput
from numpy import zeros, array, arange, linspace,float64
from scipy.integrate import simps
import convert

import aircraft_FW


class MissionResults:
    def __init__(self):
        self.time          = 0.0
        self.distance      = 0.0
        self.altitudeStart = 0.0
        self.altitudeEnd   = 0.0
        self.fuelStart     = 0.0
        self.fuelEnd       = 0.0
        self.fuelBurned    = 0.0
#        self.aircraftMassStart = 0.0
#        self.aircraftMassEnd   = 0.0
    def __repr__(self):
        out = ''
        for attr, value in self.__dict__.iteritems():
            if type(value) is float or type(value) is float64:
                out += '{:<15} = {:<12.5f}\n'.format(attr,value)
        return out

    def display(self,name=None):
        if not name==None:
            print name + '\n'+'='*len(name)
        print self.__repr__()


class Cruise(SteadyLevelFlight):
    def run_maximum_range(self,altitude,startFuel,endFuel,nSeg=10):
        fuelMass = linspace(startFuel,endFuel,nSeg)
        df = fuelMass[0]-fuelMass[1]
        totalTime = 0.0
        totalDist = 0.0
        for i,m in enumerate(fuelMass[:-1]):
            self.bi.update_fuel_mass(m)
            results = self.run_max_SAR(altitude)
            dist = results.SAR*df
            totalDist += dist
            time = dist/results.velocity
            totalTime += time
        results = MissionResults()
        results.distance      = totalDist
        results.time          = totalTime
        results.fuelStart     = startFuel
        results.fuelEnd       = endFuel
        results.fuelBurned    = startFuel - endFuel
        results.altitudeEnd   = altitude
        results.altitudeStart = altitude
        #results.aircraftMassStart = m0
        return results


    def run_maximum_endurance_time(self,altitude,startFuel,totalDuration,nSeg=10):
        time = linspace(0,totalDuration,nSeg)
        dt = time[1]-time[0]
        totalFuel = 0.0
        totalDist = 0.0
        mf = startFuel
        for i in range(nSeg-1):
            self.bi.update_fuel_mass(mf)
            results = self.run_min_fuel(altitude)
            dist = results.velocity*dt
            fuel = results.fuelFlow*dt
            totalFuel += fuel
            totalDist += dist
            mf += -fuel
        
        results = MissionResults()
        results.altitudeEnd = altitude
        results.altitudeStart = altitude
        results.distance = totalDist
        results.fuelStart = startFuel
        results.fuelBurned = startFuel - mf
        results.fuelEnd = mf
        results.time = totalDuration
        return results
    
    def run_maximum_speed(self,altitude,startFuel,endFuel):
        """ attack at sea level """
        pass
    
    def run_distance_at_speed(self,altitude,speed,distance,startFuel,nSeg=10):
        results = MissionResults()
        results.altitudeEnd = altitude
        results.altitudeStart = altitude
        results.distance = distance
        results.fuelStart = startFuel

        dist = linspace(0,distance,nSeg)
        ddist = dist[1]-dist[0]
        mf = 0.0
        dt = ddist/speed
        results.time = distance/speed
        for i in range(nSeg-1):
            self.bi.update_fuel_mass(startFuel-mf)
            self.set_results(speed,altitude)
            dmf = self.results.fuelFlow*dt
            mf += dmf
        results.fuelBurned = mf
        results.fuelEnd = startFuel-mf
        return results


class Climb(ClimbDescent):
    def run_climb(self,startAltitude,endAltitude,nSeg=10):
        altitude = linspace(startAltitude,endAltitude,nSeg)
        dh = altitude[1]-altitude[0]
        rateOfClimb = zeros(nSeg-1)
        fuelMass = zeros(nSeg-1)
        time = zeros(nSeg-1)
        distance = zeros(nSeg-1)
        totalDist = 0.0
        totalTime = 0.0
        totalFuel = 0.0
        for i,alt in enumerate(altitude[:-1]):
            data = self.run_max_climb_rate(alt)
            rateOfClimb[i] = data.climbRate
            time[i] = dh/rateOfClimb[i]
            distance[i] = data.velocity*time[i]
            fuelMass[i] = data.fuelFlow*time[i]
            totalDist += distance[i]
            totalTime += time[i]
            totalFuel += fuelMass[i]
            self.bi.update_fuel_mass_burned(fuelMass[i])
        
        results = MissionResults()
        results.altitudeEnd = endAltitude
        results.altitudeStart = startAltitude
        results.distance = totalDist
        results.time = totalTime
        results.fuelBurned = totalFuel
        return results


def run_mission_B15():
    ac = aircraft_FW.load('Baseline1')
    #ac.display()
    # -- mission inputs --
    altField = 0.0
    altCruise = 1.0e4
    altAttack = convert.ft_to_m(2000)
    distAttack = convert.nm_to_m(200.0)
    timeReserve = 1800.0 # 30min
    speedMaxAttack = convert.kt_to_msec(360)
    
    Wf0 = ac.mass.fuel.mass
    # -- assumptions --
    fuelReserveStart = 0.2*Wf0
    fuelAttackStart = 0.5*Wf0
    
    # -- calculations --
    slf = Cruise(ac)
    clm = Climb(ac)
    
    # climb 1
    climb1 = clm.run_climb(altField, altCruise)
    climb1.display('Climb 1')

    # reserve fuel
    reserve = slf.run_maximum_endurance_time(altField,fuelReserveStart,timeReserve)
    reserve.display('Reserve')
    WfReserve = reserve.fuelBurned + 0.05*Wf0
    
    # penetration-withdrawal
    slf.bi.ac.mass.display()
    penetration = slf.run_distance_at_speed(altAttack,speedMaxAttack,distAttack,fuelAttackStart)
    penetration.display('Penetration')
    slf.drop_payload()
    withdrawal = slf.run_distance_at_speed(altAttack,speedMaxAttack,distAttack,penetration.fuelEnd)
    withdrawal.display('Withdrawal')
    slf.bi.ac.mass.display()

def _run_mission15():
    ac = aircraft_FW.load('X45C')
    bi = BasicInput(ac)
    tm = ac.propulsion
    slf = Cruise(bi,tm)
    clm = Climb(bi,tm)
    
    hCruise = 10000 # cruise altitude, m
    hPenetration = convert.ft_to_m(2000.) # penetration altitude, m
    hField = 0.0
    hReserve = 0.0
    reserveTime = 30.0 #min
    attackDistance = 370400. #m -> 200nm
    
    reserveTime *= 60.
    
    climb = clm.run_climb(hField, hCruise)
    m1 = ac.mass.fuel.mass
    distance =  slf.run_maximum_range(hCruise,m1,climb[2])
    print m1
    print m1-1000.
    print distance
    reserve = slf.run_min_fuel(hReserve)
    reserveFuel = reserve.fuelFlow*reserveTime


def run_mission11():
    ac = aircraft_FW.load('Baseline1')
    #ac.display()
    #ac.mass.display()
    #bi = BasicInput(ac)
    #tm = ac.propulsion
    slf = Cruise(ac)
    clm = Climb(ac)
    
    hCruise = 10000 # cruise altitude, m
    hField = 0.0
    hReserve = 0.0
    reserveTime = 30.0 #min

    reserveTime *= 60.
    
    distClm, timeClm, fuelClm = clm.run_climb(hField, hCruise)
    reserve = slf.run_min_fuel(hReserve)
    reserveFuel = reserve.fuelFlow*reserveTime
    m1 = ac.mass.fuel.mass
    distCrs, timeCrs =  slf.run_maximum_range(hCruise,m1,reserveFuel)
    print m1
    print m1-1000.
    print distCrs/1e3
    print timeCrs/3600
    reserve = slf.run_min_fuel(hReserve)
    reserveFuel = reserve.fuelFlow*reserveTime


if __name__=="__main__":
    run_mission_B15()