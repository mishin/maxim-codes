# -*- coding: utf-8 -*-
"""
Created on Fri Sep 05 14:08:07 2014

@author: Maxim
"""

from cfd_wing import *
from misc_tools import read_tabulated_data_without_header
from ucav_design_1 import DesignFormulation

def create_input_files():
    # Longitudinal (half wing)
    DOE = read_tabulated_data_without_header('LHS_dvar9_sample30.txt')
    ac = DesignFormulation()
    ac.setup()
    ac.load_xls('Baseline1')
    
    Mach = ac.designGoals.cruiseMach
    altitude = ac.designGoals.cruiseAltitude
    yplus = 0.5
    wdir = 'D:\FW_LHS30'
    alpha = [0,2]
    beta = 2.0
    niter = 5000
    
    batFileLongitudinal = '%s\\batchLongitudinal.bat'%wdir
    fid = open(batFileLongitudinal,'at')
    pathFluent = r'C:\Program Files\Ansys Inc\v130\fluent\ntbin\win64\fluent.exe'
    for i,xnorm in enumerate(DOE):
        name = '%s\\case%d'%(wdir,i+1)
        print name
        ac.set_x(xnorm)
        ac.save_2d_figure('%s.png'%name)
        igsPath = '%s.igs'%name
        symCasPath = '%s_sym.cas'%name
        nonsymCasPath = '%s_half.cas'%name
        glfPath = '%s.glf'%name
        #create_catia_wing(ac,igsPath)
        #create_fw_cmesh(ac,igsPath,symCasPath,nonsymCasPath,yplus,glfPath)
        Sref = ac.wing.area
        Cref = ac.wing.MAC
        CG = ac.get_cg()
        for a in alpha:
            outPrefix = '%s\\results\\case%d_half_a%d'%(wdir,i+1,a*100)
            journalFileRel = 'case%d_half_a%d.jou'%(i+1, a*100)
            journalFile = '%s_half_a%d.jou'%(name, a*100)
            run_wing_half(nonsymCasPath, outPrefix, ac.designGoals.fc, a, 0, 
                          Sref, Cref, CG,journalFile,niter)
            fid.write('\"%s\" 3ddp -t7 -i %s -wait\n'%(pathFluent, journalFileRel))
    fid.close()
            

if __name__=="__main__":
    create_input_files()