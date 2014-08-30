# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 15:21:54 2014

@author: Maxim
"""
import aircraft_FW
import numpy as np
import os


def create_fw_cmesh(ac=None,igsPath=None,casPathSym=None, casPathNonsym=None, yplus=0.5):
    if ac==None:
        ac = aircraft_FW.load('Baseline1')
    if igsPath==None:
        igsPath = r'E:/1. Current work/2014 - UAV performance code/CFD automation/fw3.igs'
    if casPathSym==None:
        casPathSym = 'fw_case_sym.cas'
    if casPathNonsym==None:
        casPathNonsym = 'fw_case_nonsym.cas'

    fid = open('./templates/fw_automation.glf')
    lines = fid.readlines()
    fid.close()
    #print ac.get_cg()

    nWallChord = 65
    nWallSeg1  = 30
    nWallSeg2  = 31
    nWallTip   = 8
    nFFx = 51
    nFFy = 100
    nFFz = 55
    
    #yplus = 0.5
    dsLE = 1e-3

    dsFFy = 1e-3
    dsTopAft = 1e-2
    dsTopFront = 1e-2
    dsZ = 1e-3
    
    # igs input file
    
    lines[14] = '  $_TMP(mode_1) initialize -type Automatic {%s}\n'%igsPath
    # wing tip
    tc = ac.wing.airfoils[-1].thickness
    tc *= ac.wing.chords[-1] / 2.
    sweep = ac.wing.segSweepLErad[-1] + np.radians(5)
    pt1offset = np.array([tc*np.tan(sweep), 0, tc])
    pt2offset = np.array([0, 0, tc])
    
    lines[40] = '  $_TMP(PW_2) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 0]] {%.6f %.6f %.6f}]\n'%(pt1offset[0], pt1offset[1],pt1offset[2])
    lines[53] = '  $_TMP(PW_3) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 1]] {%.6f %.6f %.6f}]\n'%(pt2offset[0], pt2offset[1],pt2offset[2])
    
    # ff aft
    span = ac.wing.span + tc
    lines[186] = '  $_TMP(PW_13) addPoint [pwu::Vector3 add [$_CN(11) getPosition -arc 1] {0 0 %.6f}]\n'%span
    lines[198] = '  $_TMP(PW_14) addPoint [pwu::Vector3 add [$_CN(13) getPosition -arc 1] {0 0 %.6f}]\n'%span
    lines[211] = '  $_TMP(PW_15) addPoint [pwu::Vector3 add [$_CN(12) getPosition -arc 1] {0 0 %.6f}]\n'%span

    # wall spacing BL
    ds1 = ac.designGoals.fc.get_wall_spacing(yplus)
    cr = ac.wing.chords[0]
    ct = ac.wing.chords[-1]
    lines[453] = '  $_TMP(PW_38) setBeginSpacing %.6e\n'%(ds1*cr)
    lines[456] = '  $_TMP(PW_39) setBeginSpacing %.6e\n'%(ds1*cr)
    lines[464] = '  $_TMP(PW_40) setBeginSpacing %.6e\n'%(ds1*ct)
    lines[467] = '  $_TMP(PW_41) setBeginSpacing %.6e\n'%(ds1*ct)

    # leading/trailing edge spacing
    c2 = ac.wing.chords[1]
    dsLEr = dsLE*cr
    dsLE2 = dsLE*c2
    dsLEt = dsLE*ct
    
    lines[551] = '  $_TMP(PW_59) setBeginSpacing %.10f\n'%(dsLEr)
    lines[554] = '  $_TMP(PW_60) setEndSpacing %.10f\n'%(dsLEr)
    lines[557] = '  $_TMP(PW_61) setBeginSpacing %.10f\n'%(dsLEr)
    lines[560] = '  $_TMP(PW_62) setEndSpacing %.10f\n'%(dsLEr)
    
    lines[568] = '  $_TMP(PW_63) setBeginSpacing %.10f\n'%(dsLE2)
    lines[571] = '  $_TMP(PW_64) setEndSpacing %.10f\n'%(dsLE2)
    lines[574] = '  $_TMP(PW_65) setBeginSpacing %.10f\n'%(dsLE2)
    lines[577] = '  $_TMP(PW_66) setEndSpacing %.10f\n'%(dsLE2)
    
    lines[585] = '  $_TMP(PW_67) setBeginSpacing %.10f\n'%(dsLEt)
    lines[588] = '  $_TMP(PW_68) setEndSpacing %.10f\n'%(dsLEt)
    lines[591] = '  $_TMP(PW_69) setBeginSpacing %.10f\n'%(dsLEt)
    lines[594] = '  $_TMP(PW_70) setEndSpacing %.10f\n'%(dsLEt)
    lines[597] = '  $_TMP(PW_71) setBeginSpacing %.10f\n'%(dsLEt)
    lines[600] = '  $_TMP(PW_72) setEndSpacing %.10f\n'%(dsLEt)
    lines[608] = '  $_TMP(PW_73) setEndSpacing %.10f\n'%(dsLEt)

    # case file export
    lines[1094] = '  $_TMP(mode_1) initialize -type CAE {%s}\n'%casPathNonsym
    lines[1118] = '  $_TMP(mode_4) initialize -type CAE {%s}\n'%casPathSym
    
    glfPath = os.path.abspath('temp/new.glf')
    fid = open(glfPath,'wt')
    for line in lines:
        fid.write(line)
    fid.close()
    os.system('\"%s\"'%glfPath)
    os.remove(glfPath)


if __name__=="__main__":
    create_fw_cmesh()