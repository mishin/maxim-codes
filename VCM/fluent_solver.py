# -*- coding: utf-8 -*-
"""
Created on Wed Nov 06 19:53:14 2013

@author: Maxim
"""
from numpy import arange, radians, sin, cos, zeros
from pointwise_mesh import ScriptFile
from paths import CFD_paths
from os import system

class FluentOutput():
    def __init__(self):
        self.Mach = 0.0
        self.Re = 0.0
        self.alpha = 0.0
        self.cl = 0.0
        self.cd = 0.0
        self.cm = 0.0
        self.cp = list()
        self.LD = 0.0
    
    def __repr__(self):
        out = 'Mach = %.4f\n'
        out += 'Re   = %.2e\n'
        out += 'cl\tcd\tcm\n'
        out += '%.2e\t%.2e\t%.2e\n'%(self.cl,self.cd, self.cm)
        return out

class FluentAirfoil():
    def __init__(self):
        self.residuals = {'continuity':1e-3,'xvelocity':1e-3,'yvelocity':1e-3,
        'energy':1e-6,'nut':1e-3,'k':1e-3,'omega':1e-3,'epsilon':1e-3}
        self.relaxationFactor = {'a':0.1}
        self.turbulenceInput = {'SA':'no\nno\nyes\nno\n10\n',
                                'ke-realizable':'no\nno\nyes\n10\n10\n'}
        self.residualsInput = {'SA':['continuity','xvelocity','yvelocity','energy','nut'],
                               'ke-realizable':['continuity','xvelocity','yvelocity','energy','k','epsilon']}
        self.turbulenceName = {'SA':'spalart-allmaras','ke-realizable':'ke-realizable'}
        self.iterMax = 5000
        self.momentAxis = [0.25,0]
        self.turbulence = 'SA' #'ke-realizable'
        self.paths = CFD_paths()
        self.result = FluentOutput()
    
    def _create_journal_file(self,alpha,flightConditions,caseFilePath=None,turbulenceModel='SA',Cp=False,
                             journalPath=None,outputDirectory=None):
        ffbc = 'bc-3-5'
        wallbc = 'bc-2-4'
        if outputDirectory==None:
            outputDirectory = self.paths.tmpdir
        if journalPath==None:
            journalPath = self.paths.file_jou
        if caseFilePath==None:
            caseFilePath = self.paths.file_cas
        self.paths.set_name_alpha(alpha)
        freestream = [cos(radians(alpha)), sin(radians(alpha))]
        script = open(journalPath,'wt')
        script.write('/file/read\n')
        script.write('%s\n'%caseFilePath)
        script.write('/define/operating-conditions/operating-pressure\n')
        script.write('0\n')
        script.write('/define/models/viscous/%s\nyes\n'%self.turbulenceName[turbulenceModel])
        if flightConditions.Mach>=0.7:
            script.write('/define/models/solver/density-based-implicit\nyes\n')
        script.write('/define/materials/change-create\n')
        script.write('air\n\nyes\n')
        script.write('ideal-gas\n')
        script.write('no\nno\nyes\n')
        script.write('sutherland\n')
        script.write('three-coefficient-method\n')
        script.write('1.716e-05\n273.11\n110.56\n')
        script.write('no\nno\nno\n')
        script.write('/define/boundary-conditions/pressure-far-field/\n')
        script.write('%s\n'%ffbc)
        script.write('no\n')
        script.write('%.2f\n'%flightConditions.atmosphere.pressure)
        script.write('no\n%.6f\n'%flightConditions.Mach)
        script.write('no\n%.6f\n'%flightConditions.atmosphere.temperature)
        script.write('no\n')
        script.write('%.10f\n'%freestream[0])
        script.write('no\n')
        script.write('%.10f\n'%freestream[1])
        script.write('%s'%self.turbulenceInput[turbulenceModel])
        script.write('/solve/monitors/residual/convergence-criteria\n')
        for resid in self.residualsInput[turbulenceModel]:
            script.write('%.4e\n'%self.residuals[resid])
        script.write('/solve/monitors/force/drag-coefficient\n')
        script.write('yes\n%s\n\nno\nyes\n'%wallbc)
        script.write('\"%s\"\nno\nno\n'%self.paths.file_cd_hist)
        script.write('%.10f\n'%freestream[0])
        script.write('%.10f\n'%freestream[1])
        script.write('/solve/monitors/force/lift-coefficient\n')
        script.write('yes\n%s\n\nno\nyes\n'%wallbc)
        script.write('\"%s\"\nno\nno\n'%self.paths.file_cl_hist)
        script.write('%.10f\n'%(-freestream[1]))
        script.write('%.10f\n'%freestream[0])
        script.write('/solve/monitors/force/moment-coefficient\n')
        script.write('yes\n%s\n\nno\nyes\n'%wallbc)
        script.write('\"%s\"\nno\nno\n'%self.paths.file_cm_hist)
        script.write('%.4f\n'%self.momentAxis[0])
        script.write('%.4f\n'%self.momentAxis[1])
        script.write('0\n0\n1\n')
        script.write('/report/reference-values/compute/pressure-far-field\n%s\n'%ffbc)
        script.write('/solve/initialize/compute-defaults/pressure-far-field\n%s\n'%ffbc)
        script.write('/solve/iterate\n')
        script.write('%d\n'%self.iterMax)
        script.write('\nexit\nok\n')
        script.close()
    
    def run_at_aoa(self,alpha,flightConditions,caseFilePath=None,
                   turbulenceModel='SA',Cp=False,iterMax=5000):
        """
        Run Ansys fluent airfoil analysis at single angle of attack
        
        Parameters
        ----------
        alpha : float, deg
            angle of attack
        flightCondtions : object
            flight condtions object
        caseFilePath : string
            path of the case file with mesh
        turbulenceModel : string
            two models are available now "SA" and "ke-realizable"
        Cp : bool
            output Cp distribution. Not available
        iterMax : int
            maximum number of iterations
        """
        self.result.alpha = alpha
        self.result.Mach = flightConditions.Mach
        self.result.Re = flightConditions.Re
        self.iterMax = iterMax
        self._create_journal_file(alpha,flightConditions,caseFilePath,turbulenceModel,Cp)
        self._run_fluent()
        self._collect_output()
        return self.result

    def _run_fluent(self):
        system('\"\"%s\" 2ddp -hidden -i \"%s\"\"'%(self.paths.fluent,self.paths.file_jou))
    
    def _collect_output(self,histFileDir=None,Cp=False):
        if histFileDir==None:
            histFileDir=self.paths.tmpdir
        self._read_history_files(self.paths.list_hist_files)
    
    def _read_history_files(self,listOfHistFilesPath):
        result = zeros(3)
        for i,histFilePath in enumerate(listOfHistFilesPath):
            result[i] = self._read_history_file(histFilePath)
        self.result.cl = result[0]
        self.result.cd = result[1]
        self.result.cm = result[2]
        self.result.LD = result[0]/result[1]

    def _read_history_file(self,histFilePath):
        fid = open(histFilePath,'rt')
        line = fid.readlines()[-1]
        fid.close()
        return float(line.split()[1])