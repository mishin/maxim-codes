# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 22:55:54 2013

@author: Maxim
"""
from numpy import array, arange, flipud, hstack, vstack
from pointwise_mesh import AirfoilCMesh, AirfoilOMesh
from paths import CFD_paths
import airfoil
from FlightConditions import FlightConditions
from fluent_solver import FluentAirfoil, FluentOutput

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
        self.fluent._meshtype = mesh
    
    def create_mesh(self,reverse=False):
        self.mesh.create(self.paths.file_glf,self._caseFilePath)

    def run_for_multiple_aoa(self,alphaSequence=arange(-5,20,3),
                             turbulenceModel='SA',Cp=False,iterMax=5000,densityBased=False):
        """
        Run fluent at sequence of angle of attacks
        
        Parameters
        ----------
        
        alphaSequence : array, deg
            array of angles of attack
        turbulenceModel : string
            fluent turbulence model. Two options are available: 'SA' and 
            'ke-realizable'
        Cp : bool
            This option is not available yet.
        iterMax : int
            maximum number of iterations
        """
        result = FluentOutput(len(alphaSequence))
        for i,alpha in enumerate(alphaSequence):
            result1 = self.fluent.run_at_aoa(alpha,self.fc,self._caseFilePath,turbulenceModel,Cp,iterMax,densityBased)
            result.alpha[i] = result1.alpha
            result.cl[i] = result1.cl
            result.cd[i] = result1.cd
            result.cm[i] = result1.cm
            result.LD[i] = result1.LD
        return result
    
    def run_for_single_aoa(self,alpha=0.0,turbulenceModel='SA',Cp=False,iterMax=5000,densityBased=False):
        """
        Run fluent at single angle of attack
        
        Parameters
        ----------
        
        alpha : float, deg
            angle of attack
        turbulenceModel : string
            fluent turbulence model. Two options are available: 'SA' and 
            'ke-realizable'
        Cp : bool
            This option is not available yet.
        iterMax : int
            maximum number of iterations
        """
        result = self.fluent.run_at_aoa(alpha,self.fc,self._caseFilePath,turbulenceModel,Cp,iterMax,densityBased)
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

def run_lowRe():
    af = airfoil.Airfoil()
    af.read_txt(r'D:\laptop_sync\1. Projects\RENN\airfoil design\AG24new.txt')
#    result = af.get_X_polar(0.02,0.355e6,[0,16,2.0],nIter=100)
    fc = FlightConditions(22.,0,0,0.24)
    solver = CFDsolver(af,fc,1.0,mesh='O')
    solver.fluent.residuals['energy']=1e-6
    solver.fluent.relaxationFactor['xvelocity'] = 1e-3
    solver.mesh._airfoilPts = 75
    solver.mesh._interiorPts = 75
    solver.mesh._dsTE = 1e-5
    solver.mesh._dsLE = 2e-3
    solver.mesh._growthRate = 1.2
    solver.create_mesh()
    alpha = array([-10.,-5,0,5,10,12,14,16,18])
    result = solver.run_for_multiple_aoa(alpha,turbulenceModel='SA')
    for a,cl,cd,cm in zip(result.alpha,result.cl,result.cd,result.cm):
        print a,'\t',cl,'\t',cd,'\t',cm

if __name__=="__main__":
    run_lowRe()