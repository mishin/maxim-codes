# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 22:55:54 2013

@author: Maxim
"""
from numpy import array, arange
from pointwise_mesh import AirfoilMesh
from paths import CFD_paths
import airfoil
from FlightConditions import FlightConditions
from fluent_solver import FluentAirfoil

class CFDsolver():
    def __init__(self,airfoil,flightConditions,yplus=1.0):
        self.af = airfoil
        self.fc = flightConditions
        self.paths = CFD_paths()
        self._caseFilePath = self.paths.file_cas
        self.mesh = AirfoilMesh(airfoil,dSwall=self.fc.get_wall_spacing(yplus))
        self.mesh.create(self.paths.file_glf,self._caseFilePath)
        self.fluent = FluentAirfoil()
    
    def run_for_multiple_aoa(self,alphaSequence=arange(-5,20,3),Cp=False):
        for alpha in alphaSequence:
            self.fluent.run_at_aoa(alpha,self.fc,self._caseFilePath,Cp)
    
    def run_for_single_aoa(self,alpha=0.0):
        self.fluent.run_at_aoa(alpha,self.fc,self._caseFilePath,Cp)

def run_debug1():
    af = airfoil.Airfoil()
    #af.read_txt('GA37A315.txt')
    af.naca4()
    fc = FlightConditions(250,3e3)
    fc.atmosphere.pressure
    solver = CFDsolver(af,fc)
    solver.fluent.residuals['energy']=1e-3
    solver.fluent._create_journal_file(3,fc,turbulenceModel='ke-realizable')
    

if __name__=="__main__":
    run_debug1()