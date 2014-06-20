# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 16:50:41 2014

@author: Maxim
"""

from performance import SteadyLevelFlight, ClimbDescent, BasicInput
from numpy import zeros, array, arange, linspace
from scipy.integrate import simps

class MissionResults:
    """
    stores results of mission analysis such as velocity, mass, power setting etc. 
    at different segments of the mission.
    """
    def __init__(self,nseg=10):
        self.nseg = nseg
        self.n    = 0
        self.totalMass       = zeros(nseg)
        self.fuelMass        = zeros(nseg)
        self.velocity        = zeros(nseg)
        self.density         = zeros(nseg)
        self.altitude        = zeros(nseg)
        self.fuelFlow        = zeros(nseg)
        self.thrust          = zeros(nseg)
        self.SAR_km_per_kg   = zeros(nseg)
        self.drag            = zeros(nseg)
        self.lift            = zeros(nseg)
        self.range           = 0.0
        self.endurance       = 0.0
        self.fuelBurned      = 0.0
    def add_performance_result(self,rslt,fuelMass=0,totalMass=0):
        i = self.n
        self.n += 1
        self.velocity[i]        = rslt.velocity
        self.EAS[i]             = rslt.equivalentAirspeed
        self.altitude[i]        = rslt.densityAltitude
        self.density[i]         = rslt.density
        self.fuelFlow[i]        = rslt.fuelFlow
        self.thrust[i]          = rslt.thrust
        self.power[i]           = rslt.power
        self.SAR_km_per_kg[i]   = rslt.SAR_km_per_kg
        self.drag[i]            = rslt.drag
        self.lift[i]            = rslt.lift
        self.fuelMass[i]        = fuelMass
        self.totalMass[i]       = totalMass
        if self.nseg>=2:
            self._integrate()
    def _integrate(self):
        self.range = simps(self.SAR_km_per_kg[:self.n],self.fuelMass[:self.n])
        specificEndurance = 1.0/(array(self.fuelFlow[:self.n])*3600.0)
        self.endurance = simps(specificEndurance,self.fuelMass[:self.n])
        self.fuelBurned = abs(self.fuelMass[0] - self.fuelMass[self.n])



class Cruise(SteadyLevelFlight):
    def run_maximum_range(self,altitude,startFuel,endFuel):
        pass
    
    def run_maximum_endurance(self,altitude,startFuel,endFuel):
        pass
    
    def run_maximum_speed(self,altitude,startFuel,endFuel):
        pass


class Climb(ClimbDescent):
    def run_climb(self,startAltitude,endAltitude,nSeg):
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
            
        return totalDist, totalTime, totalFuel


def run_test1():
    import aircraft_FW
    
    ac = aircraft_FW.load('X45C')
    bi = BasicInput(ac)
    clm = Climb(bi,ac.propulsion)
    print clm.run_climb(0,10000,10)


if __name__=="__main__":
    run_test1()