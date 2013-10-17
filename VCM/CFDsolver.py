# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 15:03:44 2012

@author: Maxim
"""

import numpy as ny
import os
from math import sin, cos, pi, radians

from ISAtmosphere import ISAtmosphere
import paths
from Airfoil import Airfoil

class Flight_conditions():
    def __init__(self, altitude = 0, velocity = 50):
        self.velocity = velocity
        self.ref_length = 1.0
        self.ISA = ISAtmosphere(altitude)
        self.mu = self.ISA.sutherlandVisc()
        self.Re = self.ISA.density * velocity * self.ref_length / self.mu
        self.Mach = velocity / self.ISA.soundSpeed

class Airfoil_mesh():
    def __init__(self,paths, flight_conditions, create_glf = True, run = True):
        self.Xmax = 20.
        self.Ymax = 15.
        self.XtopOffset = -20.
        self.dsLE = 1e-3
        self.dsTE = 1e-3
        self.dsFarfield = 1.
        self.dsTop = 1e-2
        self.afPt = 101
        self.Ypt = 75
        self.Xpt = 101
        self.yplus_wall = 1.
        self.flight_cond = flight_conditions
        
        if run: create_glf = True
        if create_glf: self._create_glf(paths)
        if run:
            cmd = ('\"%s\"'%paths.file_glf)
            os.system(cmd)

    def ds(self,yplus):
        fc = self.flight_cond
        Cf = 0.026 / (fc.Re**(1/7))
        tau_wall = Cf * fc.ISA.density * fc.velocity**2/2
        U_fric = (tau_wall/fc.ISA.density)**0.5
        return yplus*fc.mu/(U_fric*fc.ISA.density)

    def _create_glf(self,paths):
        inData = ('$_TMP(mode_10) initialize -type Automatic {%s}'%paths.file_igs)
        line = ('$_TMP(PW_57) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 1]] {%.0f 0 0}]'%self.Xmax)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_58) addPoint [pwu::Vector3 add [$_CN(3) getPosition -arc 1] {0 %.0f 0}]'%self.Ymax)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_59) addPoint [pwu::Vector3 add [$_CN(4) getPosition -arc 1] {%.0f 0 0}]'%self.XtopOffset)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_64) do setDimension %d'%self.Ypt)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_65) do setDimension %d'%self.Xpt)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_66) do setDimension %d'%(2*self.afPt-1))
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_67) do setDimension %d'%self.afPt)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_68) setBeginSpacing %.6f'%self.ds(self.yplus_wall))
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_69) setBeginSpacing %.6f'%self.ds(self.yplus_wall))
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_70) setEndSpacing %d'%self.dsFarfield)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_71) setEndSpacing %d'%self.dsFarfield)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_72) setBeginSpacing %d'%self.dsFarfield)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_73) setEndSpacing %d'%self.dsFarfield)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_74) setBeginSpacing %d'%self.dsFarfield)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_75) setBeginSpacing %.6f'%self.dsTE)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_76) setEndSpacing 0.001')
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_77) setEndSpacing 0.001')
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_78) setBeginSpacing %.6f'%self.dsTop)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_79) setEndSpacing %.6f'%self.dsTop)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_80) setBeginSpacing %.6f'%self.dsLE)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_81) setBeginSpacing %.6f'%self.dsLE)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_82) setEndSpacing %.6f'%self.dsTE)
        inData = ny.vstack([inData,line])
        line = ('$_TMP(PW_83) setEndSpacing %.6f'%self.dsTE)
        inData = ny.vstack([inData,line])
        line = ('pw::Application export [list $_DM(1)] {%s}'%paths.file_cas)
        inData = ny.vstack([inData,line])

        lineNumber = ny.array([[11,36,49,63,116,123,131,138,176,179,187,190,193,196,199,207,210,213,221,224,232,235,243,246,273]]).T
        inData = ny.hstack([lineNumber,inData])
        template = open(paths.template_pw,'rt')
        templateText = template.readlines()
        template.close()
        for line in inData: templateText[int(line[0])] = (line[1]+'\n')
        jouFile = open(paths.file_glf,'wt')
        for line in templateText: jouFile.write(line)
        jouFile.close()

class Solver():
    def __init__(self, paths, flight_cond):
        self.flight = flight_cond
        self.operating_pressure = 0
        self.turb_model = 'ke-realizable' #spalart-allmaras, kw-sst
        self.gauge_pressure = flight_cond.ISA.pressure
        self.resid_cont = 0.001
        self.resid_Xvel = 0.001
        self.resid_Yvel = 0.001
        self.resid_energy = 1e-6
        self.resid_nut = 0.001
        self.resid_k = 0.001
        self.resid_omega = 0.001
        self.resid_epsilon = 0.001
        self.momentAxis = ([0.25,0])
        self.maxIter = 5000
        self.paths = paths

        self.alpha = []
        self.cl = []
        self.cd = []
        self.cm = []
        
    def run_fluent(self, alpha):
        
        self.paths.set_name_alpha(alpha)
        self.alpha = ny.append(self.alpha, alpha)
        alpha = radians(alpha)
        drag_vect = ny.array([cos(alpha), sin(alpha)])
        lift_vect = ny.array([cos(alpha + pi/2), sin(alpha + pi/2)])
        
        residuals = '%.6f\n%.6f\n%.6f\n%.6f\n'%(self.resid_cont,self.resid_Xvel,self.resid_Yvel,self.resid_energy)
        if self.turb_model == 'ke-realizable':
            residuals += '%.6f\n%.6f'%(self.resid_k,self.resid_epsilon)
            #turb_specification = 'yes\nno\n1\nno\n1'
            turb_specification = 'no\nno\nyes\n10\n0.001' #turb visc ratio def 10(0.001)
        elif self.turb_model == 'spalart-allmaras':
            residuals += '%.6f\n'%self.resid_nut
            turb_specification = 'no\nno\nyes\nno\n10'
        elif self.turb_model == 'kw-sst': #do not use, need to debug
            residuals += '%.6f\n%.6f'%(self.resid_k, self.resid_omega)
            turb_specification = 'yes\nno\n1\nno\n1'
        
        input_data = ('\"%s\"'%self.paths.file_cas)
        line = ('%.0f'%self.operating_pressure)
        input_data = ny.vstack([input_data,line])
        line = ('/define/models/viscous/%s'%self.turb_model)
        input_data = ny.vstack([input_data,line])
        line = ('%.0f'%self.gauge_pressure)
        input_data = ny.vstack([input_data,line])
        line = ('%.4f'%self.flight.Mach)
        input_data = ny.vstack([input_data,line])
        line = ('%.1f'%self.flight.ISA.temperature)
        input_data = ny.vstack([input_data,line])
        line = ('%.6f'%drag_vect[0])
        input_data = ny.vstack([input_data,line])
        line = ('%.6f'%drag_vect[1])
        input_data = ny.vstack([input_data,line])
        input_data = ny.vstack([input_data,turb_specification])
        input_data = ny.vstack([input_data, residuals])
        line = ('\"%s\"'%self.paths.file_cd_hist)
        input_data = ny.vstack([input_data,line])
        line = ('%.6f'%drag_vect[0])
        input_data = ny.vstack([input_data,line])
        line = ('%.6f'%drag_vect[1])
        input_data = ny.vstack([input_data,line])
        line = ('\"%s\"'%self.paths.file_cl_hist)
        input_data = ny.vstack([input_data,line])
        line = ('%.6f'%lift_vect[0])
        input_data = ny.vstack([input_data,line])
        line = ('%.6f'%lift_vect[1])
        input_data = ny.vstack([input_data,line])
        line = ('\"%s\"'%self.paths.file_cm_hist)
        input_data = ny.vstack([input_data,line])
        line = ('%.4f'%self.momentAxis[0])
        input_data = ny.vstack([input_data,line])
        line = ('%.4f'%self.momentAxis[1])
        input_data = ny.vstack([input_data,line])
        line = ('%d'%self.maxIter)
        input_data = ny.vstack([input_data,line])
                
        lineNumber = ny.array([[1,3,4,25,27,29,31,33,34,36,43,46,47,54,57,58,65,68,69,78]]).T
        input_data = ny.hstack([lineNumber,input_data])
        
        template = open(self.paths.template_fl,'rt')
        templateText = template.readlines()
        template.close()
        for line in input_data: templateText[int(line[0])] = (line[1]+'\n')

        jouFile = open(self.paths.file_jou,'wt')
        for line in templateText: jouFile.write(line)
        jouFile.close()
        
        cmd = ('\"\"%s\" 2ddp -hidden -i \"%s\"\"'%(self.paths.fluent,self.paths.file_jou))
        #cmd = ('\"\"%s\" 2ddp -i \"%s\"\"'%(self.paths.fluent,self.paths.file_jou))
        print cmd
        os.system(cmd)
        
        hist_files = [self.paths.file_cl_hist, self.paths.file_cd_hist, self.paths.file_cm_hist]
        coef = []
        for hist_file in hist_files:
            histFile = open(hist_file,'rt')
            lines = histFile.readlines()
            histFile.close()
            line = lines[-1]
            coef.append( float(line.split()[1]) )
        self.cl = ny.append(self.cl, coef[0])
        self.cd = ny.append(self.cd, coef[1])
        self.cm = ny.append(self.cm, coef[2])

def testFcn():
    path = paths.CFD_paths('GA_opt')
    cruise = Flight_conditions(1500, 48.92)
    af = Airfoil()
    af.CST_airfoil(ny.array([0.178551,0.273254,0.268906,0.226346]),ny.array([-0.178551,-0.101338,-0.255260,-0.043527]))
    af.write_airfoil_txt('optAirfoil.dat')
    af.create_af_CAT(save = path.file_igs)
    Airfoil_mesh(path,cruise)
    Airfoil_mesh.yplus_wall = 1.0
    fluent = Solver(path,cruise)
    fluent.turb_model = 'ke-realizable'
    for alpha in range(10, 21, 2):
        fluent.run_fluent(alpha)
    af.polar.Re = cruise.Re
    af.polar.M = cruise.Mach
    af.polar.alpha = fluent.alpha
    af.polar.CL = fluent.cl
    af.polar.CD = fluent.cd
    af.polar.CM = fluent.cm
    #af.calc_Jpolar(0.1463,2999951,ny.array([-2,20,2.]))
    af.write_polar_txt(r'D:\3. Projects\afOpt\experimental data\GAopt.pol')

#def batch_calculation():
    designs = open(r'D:\3. Projects\afOpt\results\CFD_20120719\Designs.txt','rt')
    Au = ny.zeros([4,1])
    Al = ny.zeros([4,1])
    lines = designs.readlines()
    designs.close()
    landing = Flight_conditions(0, 32.5)
    path = paths.CFD_paths()
    results = r'D:\3. Projects\afOpt\results\CFD_20120719'
    for N,line in enumerate(lines):
        path.set_name(('DoE_%d'%N))
        segLine = line.split()
        Au[0] = float(segLine[0])
        Au[1] = float(segLine[1])
        Au[2] = float(segLine[2])
        Au[3] = float(segLine[3])
        Al[0] = -float(segLine[0])
        Al[1] = float(segLine[4])
        Al[2] = float(segLine[5])
        Al[3] = float(segLine[6])

        af = Airfoil()
        af.name = ('DoE %d'%N)
        af.CST_airfoil(Au,Al)
        af.create_af_CAT(save = path.file_igs)
        Airfoil_mesh(path,landing)
        Airfoil_mesh.yplus_wall = 1.0
        fluent = Solver(path,landing)
        fluent.turb_model = 'ke-realizable'
        alphaSeq = ny.array([10, 12, 14, 16, 18, 20])
        for alpha in alphaSeq:
            fluent.run_fluent(alpha)
        af.polar.Re = landing.Re
        af.polar.M = landing.Mach
        af.polar.alpha = fluent.alpha
        af.polar.CL = fluent.cl
        af.polar.CD = fluent.cd
        af.polar.CM = fluent.cm
        af.write_airfoil_txt(('%s\\airfoil_%d_coord.dat'%(results,N)))
        af.write_polar_txt(('%s\\airfoil_%d.pol'%(results,N)))

def calc_airfoil(name,polarPath):
#    Au = ny.array([0.178551,0.273254,0.268906,0.226346])
#    Al = ny.array([-0.178551,-0.101338,-0.255260,-0.043527])
    Au = ny.array([0.1920,0.2887,0.2844,0.2304])
    Al = ny.array([-0.1920,-0.0859,-0.2398,-0.0255])
    path = paths.CFD_paths()
    landing = Flight_conditions(1500.0,60.0)
    af = Airfoil()
    #af.read_airfoil_txt(name)
    af.CST_airfoil(Au,Al)

    af.create_af_CAT(save = path.file_igs)
    Airfoil_mesh(path,landing)
    Airfoil_mesh.yplus_wall = 1.0
    fluent = Solver(path,landing)
    fluent.turb_model = 'ke-realizable'
    #alphaSeq = ny.linspace(0,20,21)
    alphaSeq = ny.linspace(0,10,5)
    for alpha in alphaSeq:
        fluent.run_fluent(alpha)
    af.polar.Re = landing.Re
    af.polar.M = landing.Mach
    af.polar.alpha = fluent.alpha
    af.polar.CL = fluent.cl
    af.polar.CD = fluent.cd
    af.polar.CM = fluent.cm

    #af.write_polar_txt(r'D:\3. Projects\afOpt\results\af_GA20120716\landingCFDpolar_01_.pol')
    af.write_polar_txt(polarPath)

#calc_airfoil()
#batch_calculation()
#testFcn()

def calc_Re():
    v_cruise = 217 * 5/18
    v_stall = 19
    cruise = Flight_conditions(1500,v_cruise)
    landing = Flight_conditions(0, 1.2*v_stall)
    cruise.ISA.gas = 'ideal-gas'
    print v_cruise, landing.Re, landing.Mach

def run_batch_txt():
    wdir = r'D:\light aircraft\airfoil selection\CFD_2013'
    af = ['GA37A315.dat','GA37A315mod_reworked.dat','NACA747A315.dat']
    for a in af:
        path = wdir+'\\'+a
        print path
        calc_airfoil(path,path+'polarc81')

def batch_doe_txt():
    doePath = 'DoE_desings_20130603.txt'
    fid = open(doePath,'rt')
    for i,line in enumerate(fid):
        seg = line.split()
        lineData = ny.zeros([len(seg)])
        for j,val in enumerate(seg):
            lineData[j] = float(val)
        if i==0:
            data = lineData
        else:
            data = ny.vstack([data,lineData])
    fid.close()
    
if __name__=="__main__":
    calc_airfoil('NACA4415.txt','NACA4415polar.txt')
#calc_Re()