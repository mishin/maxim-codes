# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 22:55:54 2013

@author: Maxim
"""
from numpy import array, arange
from pointwise_mesh import AirfoilCMesh, AirfoilOMesh
from paths import CFD_paths
import airfoil
from FlightConditions import FlightConditions
from fluent_solver import FluentAirfoil

class CFDsolver():
    def __init__(self,airfoil,flightConditions,yplus=1.0,mesh='C'):
        self.af = airfoil
        self.fc = flightConditions
        self.paths = CFD_paths()
        self._caseFilePath = self.paths.file_cas
        if mesh=='C':
            self.mesh = AirfoilCMesh(airfoil,dsWall=self.fc.get_wall_spacing(yplus))
        elif mesh=='O':
            self.mesh = AirfoilOMesh(airfoil,dsWall=self.fc.get_wall_spacing(yplus))
        self.fluent = FluentAirfoil()
    
    def create_mesh(self):
        self.mesh.create(self.paths.file_glf,self._caseFilePath)

    def run_for_multiple_aoa(self,alphaSequence=arange(-5,20,3),Cp=False):
        for alpha in alphaSequence:
            self.fluent.run_at_aoa(alpha,self.fc,self._caseFilePath,Cp)
    
    def run_for_single_aoa(self,alpha=0.0,turbulenceModel='SA',Cp=False,iterMax=5000):
        result = self.fluent.run_at_aoa(alpha,self.fc,self._caseFilePath,turbulenceModel,Cp,iterMax)
        return result

def run_debug1():
    af = airfoil.Airfoil()
    #af.read_txt('GA37A315.txt')
    Au = array([0.119087477, 0.160950359,0.203634413,0.192468212])
    Al = array([-0.119087477, -0.200580639, -0.126010045, 0.107256400e-18])
    af.create_CST(Au,Al)
    fc = FlightConditions(0.0,9e3)
    V = fc.atmosphere.soundSpeed * 0.78
    fc = FlightConditions(V,9e3)
    fc.atmosphere.pressure
    solver = CFDsolver(af,fc,100)
    solver.fluent.residuals['energy']=1e-4
    solver.fluent.relaxationFactor['xvelocity'] = 1e-3
    solver.mesh._airfoilPts = 50
    solver.mesh._yfarfieldPts = 33
    solver.mesh._xfarfieldPts = 50
    solver.mesh._frontFarfieldSpacing = 0.1
    solver.mesh._farfieldSpacing = 1.0
    solver.create_mesh()
    lowFidelity = solver.run_for_single_aoa(0.0,iterMax=10000,turbulenceModel='ke-realizable')
    solver2 = CFDsolver(af,fc,10.0)
    solver2.mesh._airfoilPts = 100
    solver2.mesh._xfarfieldPts = 100
    solver2.mesh._yfarfieldPts = 100
    solver2.mesh._farfieldSpacing = 1.0
    solver2.mesh._frontFarfieldSpacing = 0.1
    solver2.create_mesh()
    highFidelity = solver2.run_for_single_aoa(0.0,iterMax=10000,turbulenceModel='ke-realizable')
    print lowFidelity
    print highFidelity

def run_o_mesh():
    array([0.12071848, 0.13529101, 0.1749148, 0.14246778,-0.15070396,-0.17601022,-0.04001307])
    af = airfoil.Airfoil()
    Au = array([0.12071848, 0.13529101, 0.1749148, 0.14246778])
    Al = array([-0.12071848, -0.15070396,-0.17601022,-0.04001307])
    af.create_CST(Au,Al)
    fc = FlightConditions(0.0,9e3)
    V = fc.atmosphere.soundSpeed * 0.73
    fc = FlightConditions(V,9e3)
    fc.atmosphere.pressure
    solver = CFDsolver(af,fc,50,mesh='O')
    solver.fluent.residuals['energy']=1e-6
    solver.fluent.relaxationFactor['xvelocity'] = 1e-3
    solver.mesh._airfoilPts = 50
    solver.mesh._interiorPts = 45
    solver.mesh._dsTE = 5e-4
    solver.mesh._dsLE = 2e-3
    solver.mesh._growthRate = 1.3
    solver.create_mesh()
    lowFidelity = solver.run_for_single_aoa(2.0,iterMax=10000,turbulenceModel='ke-realizable')

if __name__=="__main__":
    run_o_mesh()