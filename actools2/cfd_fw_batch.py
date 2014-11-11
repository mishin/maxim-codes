# -*- coding: utf-8 -*-
"""
Created on Fri Sep 05 14:08:07 2014

@author: Maxim
"""

from cfd_wing import *
from misc_tools import read_tabulated_data_without_header


def create_input_files():
    from ucav_design_1 import DesignFormulation
    # Longitudinal (half wing)
    doePath = r'E:\1. Current work\2014 - UAV performance code\Results\norm_results.txt'
    #DOE = read_tabulated_data_without_header('LHS_dvar9_sample30.txt')
    DOE = read_tabulated_data_without_header(doePath)
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    startCase = 30
#    Mach = ac.designGoals.cruiseMach
#    altitude = ac.designGoals.cruiseAltitude
    yplus = .5
    wdir = 'D:\sync_Maxim'
    alpha = [0,2]
#    beta = 2.0
    niter = 5000
    
    batFileLongitudinal = '%s\\batchLatDirectional.bat'%wdir
    fid = open(batFileLongitudinal,'at')
    pathFluent = r'C:\Program Files\Ansys Inc\v130\fluent\ntbin\win64\fluent.exe'
    for i,xnorm in enumerate(DOE):
        if i>5:
            name = '%s\\case%d'%(wdir,i+startCase)
            caseName = 'case%d'%(i+startCase)
            print name
            ac.set_x(xnorm)
            #FIXME: this part adds section for payload on top
            #ac._adjust_root_airfoil()
            ac.save_2d_figure('%s.png'%name)
            igsPath = '%s.igs'%name
            symCasPath = '%s_sym.cas'%name
            nonsymCasPath = '%s_half.cas'%name
            glfPath = '%s.glf'%name
            #create_catia_wing(ac,igsPath)
            create_fw_cmesh(ac,igsPath,symCasPath,nonsymCasPath,yplus,glfPath)
            Sref = ac.wing.area/2.0
            Cref = ac.wing.MAC
            #Lref = ac.wing.span
            ac._restore_root_airfoil()
            CG = ac.get_cg()
            #a = 2.0
    #        outPrefix = '%s\\results\\case%d_sym_a%d_b%d'%(wdir,i+1,a*100,beta*100)
    #        journalFileRel = 'case%d_sym_a%d_b%d.jou'%(i+1, a*100,beta*100)
    #        journalFile = '%s_sym_a%d_b%d.jou'%(name, a*100, beta*100)
            #run_wing_half(symCasPath, outPrefix, ac.designGoals.fc, a, beta, Sref, Lref, CG, journalFile, niter)
            for a in alpha:
                outPrefix = '%s\\results\\%s_half_a%d'%(wdir,caseName,a*100)
                journalFileRel = '%s_half_a%d.jou'%(name, a*100)
                journalFile = '%s_half_a%d.jou'%(name, a*100)
                run_wing_half(nonsymCasPath, outPrefix, ac.designGoals.fc, a, 0, 
                              Sref, Cref, CG,journalFile,niter)
            fid.write('\"%s\" 3ddp -t8 -i %s -wait\n'%(pathFluent, journalFileRel))
    fid.close()


def create_input_files2():
    import aircraft_FW
    ac = aircraft_FW.load('aeroValidation',False)
    name = 'aeroValidation'
    CG = [0.2,0,0]

    yplus = .5
    wdir = 'D:\sync_Maxim'
    alpha = [0,2,4,8,12]

    niter = 5000
    
    batFile = 'run_cfd_batch2.bat'
    pathFluent = r'C:\Program Files\Ansys Inc\v130\fluent\ntbin\win64\fluent.exe'
    
    igsPath = '%s.igs'%name
    symCasPath = '%s_sym.cas'%name
    nonsymCasPath = '%s_half.cas'%name
    glfPath = '%s.glf'%name
    
    create_fw_cmesh(ac,igsPath,symCasPath,nonsymCasPath,yplus,glfPath)
    Sref = ac.wing.area/2.0
    Cref = ac.wing.MAC
    
    fid = open(batFile,'wt')
    for a in alpha:
        outPrefix = '%s\\results\\%s_half_a%d'%(wdir,name,a*100)
        journalFileRel = '%s_half_a%d.jou'%(name, a*100)
        journalFile = '%s_half_a%d.jou'%(name, a*100)
        run_wing_half(nonsymCasPath, outPrefix, ac.designGoals.fc, a, 0, 
                      Sref, Cref, CG, journalFile,niter)
        fid.write('\"%s\" 3ddp -t8 -i %s -wait\n'%(pathFluent, journalFileRel))
    fid.close()



if __name__=="__main__":
    create_input_files2()