# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 11:28:44 2014

@author: Maxim
"""

from ucav_design_1 import DesignFormulation
from catia_wing import create_catia_wing
from pointwise_fw import create_fw_cmesh
from fluent_fw import run_wing
import paths
import os
from numpy import zeros
from flight_conditions import FlightConditions

class CFDresults:
    def __init__(self):
        self.Mach = list()
        self.alpha = list()
        self.beta = list()
        self.CL = list()
        self.CD = list()
        self.CY = list()
        self.Cm = list()
        self.Cn = list()
        self.Cl = list()
        # ref values
        self.Sref = 0.0
        self.Cref = 0.0
        self.bref = 0.0
        self.CG = zeros(3)
    
    def write_txt(self,path):
        self.path = os.path.abspath(path)
        fid = open(self.path,'wt')
        fid.write('Sref\t%.8f\n'%self.Sref)
        fid.write('Cref\t%.8f\n'%self.Cref)
        fid.write('bref\t%.8f\n'%self.bref)
        fid.write('CG\t%.8f\t%.8f\t%.8f\n'%(self.CG[0], self.CG[1], self.CG[2]))
        for Mach,alpha,beta,CL,CD,CY,Cm,Cn,Cl in zip(self.Mach, self.alpha, self.beta, self.CL, self.CD, self.CY, self.Cm, self.Cn, self.Cl):
            fid.write('%.8f\t'%Mach)
            fid.write('%.8f\t'%alpha)
            fid.write('%.8f\t'%beta)
            fid.write('%.8e\t'%CL)
            fid.write('%.8e\t'%CD)
            fid.write('%.8e\t'%CY)
            fid.write('%.8e\t'%Cm)
            fid.write('%.8e\t'%Cn)
            fid.write('%.8e\n'%Cl)
        fid.close()
    
    def append_txt(self,path=None,nval=1):
        fid = open(self.path, 'at')
        n = -nval
        for Mach,alpha,beta,CL,CD,CY,Cm,Cn,Cl in zip(self.Mach[n:], self.alpha[n:], self.beta[n:], self.CL[n:], self.CD[n:], self.CY[n:], self.Cm[n:], self.Cn[n:], self.Cl[n:]):
            fid.write('%.8f\t'%Mach)
            fid.write('%.8f\t'%alpha)
            fid.write('%.8f\t'%beta)
            fid.write('%.8e\t'%CL)
            fid.write('%.8e\t'%CD)
            fid.write('%.8e\t'%CY)
            fid.write('%.8e\t'%Cm)
            fid.write('%.8e\t'%Cn)
            fid.write('%.8e\n'%Cl)
        fid.close()
    
    def get_summary(self):
        """ calculate Cma, CLa, Cnb, Clb, k, CD0, CL0 """
        pass


def run_cfd_wing_analysis():
    Mach = 0.7
    altitude = 1e4 # m

    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    
    paths.myPaths.set_file_prefix('fw_debug')
    igsPath = paths.myPaths.get_tmp_file('igs')
    symCasPath = paths.myPaths.get_tmp_file('cas','_sym')
    nonsymCasPath = paths.myPaths.get_tmp_file('cas','_nonsym')
    
    print os.path.dirname(nonsymCasPath)
    print os.path.basename(nonsymCasPath)

    create_catia_wing(ac,igsPath)
    create_fw_cmesh(ac,igsPath,symCasPath, nonsymCasPath)
    os.remove(igsPath)
    
    # under development
    fc = FlightConditions(Mach, altitude)

    run_wing(nonsymCasPath, fc, True)

if __name__=="__main__":
    run_cfd_wing_analysis()
    