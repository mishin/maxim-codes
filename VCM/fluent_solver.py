# -*- coding: utf-8 -*-
"""
Created on Wed Nov 06 19:53:14 2013

@author: Maxim
"""
from numpy import arange, radians, sin, cos
from pointwise_mesh import ScriptFile
from paths import CFD_paths

class FluentAirfoil():
    def __init__(self):
        self.residuals = {'continuity':1e-3,'xvelocity':1e-3,'yvelocity':1e-3,
        'energy':1e-6,'nut':1e-3,'k':1e-3,'omega':1e-3,'epsilon':1e-3}
        self.turbulenceInput = {'SA':'no\nno\nyes\nno\n10\n',
                                'ke-realizable':'no\nno\nyes\n10\n0.001\n'}
        self.residualsInput = {'SA':['continuity','xvelocity','yvelocity','energy','nut'],
                               'ke-realizable':['continuity','xvelocity','yvelocity','energy','k','epsilon']}
        self.turbulenceName = {'SA':'spalart-allmaras','ke-realizable':'ke-realizable'}
        self.iterMax = 5000
        self.momentAxis = [0.25,0]
        self.turbulence = 'SA' #'ke-realizable'
        self.paths = CFD_paths()
    
    def _create_journal_file(self,alpha,flightConditions,caseFilePath=None,turbulenceModel='SA',Cp=False,
                             journalPath=None,outputDirectory=None):
        #TODO: remove template since fluent TUI commands are short
        if outputDirectory==None:
            outputDirectory = self.paths.tmpdir
        if journalPath==None:
            journalPath = self.paths.file_jou
        if caseFilePath==None:
            caseFilePath = self.paths.file_cas
        self.paths.set_name_alpha(alpha)
        freestream = [cos(radians(alpha)), sin(radians(alpha))]
        script = ScriptFile(self.paths.template_fl,journalPath)
        script.write_template_lines([0])
        script.write('%s\n'%caseFilePath)
        script.write_template_lines(arange(2,4))
        script.write('/define/models/viscous/%s\n'%self.turbulenceName[turbulenceModel])
        script.write_template_lines([5])
        if flightConditions.Mach>=0.7:
            script.write('/define/models/solver/density-based-implicit\nyes\n')
        script.write_template_lines(arange(6,25))
        script.write('%.2f\n'%flightConditions.atmosphere.pressure)
        script.write('no\n%.6f\n'%flightConditions.Mach)
        script.write('no\n%.6f\n'%flightConditions.atmosphere.temperature)
        script.write('no\n')
        script.write('%.10f\n'%freestream[0])
        script.write('no\n')
        script.write('%.10f\n'%freestream[1])
        script.write('%s'%self.turbulenceInput[turbulenceModel])
        script.write_template_lines([35])
        for resid in self.residualsInput[turbulenceModel]:
            script.write('%.4e\n'%self.residuals[resid])
        script.write_template_lines(arange(37,43))
        script.write('\"%s\"\nno\nno\n'%self.paths.file_cd_hist)
        script.write('%.10f\n'%freestream[0])
        script.write('%.10f\n'%freestream[1])
        script.write_template_lines(arange(48,54))
        script.write('\"%s\"\nno\nno\n'%self.paths.file_cl_hist)
        script.write('%.10f\n'%freestream[1])
        script.write('%.10f\n'%(-freestream[0]))
        script.write_template_lines(arange(59,65))
        script.write('\"%s\"\nno\nno\n'%self.paths.file_cm_hist)
        script.write('%.4f\n'%self.momentAxis[0])
        script.write('%.4f\n'%self.momentAxis[1])
        script.write_template_lines(arange(70,80))
        script.write('%d\n'%self.iterMax)
        script.write('\nexit\nok\n')
        script.close()
    
    def run_at_aoa(self,alpha,flightCondtions,caseFilePath=None,Cp=False):
        self._create_journal_file(alpha,flightConditions,caseFilePath,Cp)
        self._run_fluent()
        self._collect_output(histFilePath)
    
    def _run_fluent(self):
        pass
    
    def _collect_output(self,Cp=False):
        collect_coefficients()
        if Cp:
            collect_Cp()