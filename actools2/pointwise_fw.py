# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 15:21:54 2014

@author: Maxim
"""
import aircraft_FW
import numpy as np
import os
from geometry import rotate_2d


def create_wing_db(ac=None, glfPath=None):
    if ac==None:
        ac = aircraft_FW.load('Baseline1')
    if glfPath==None:
        glfPath='tmp_glf1.glf'
    def get_section_data(wing, idx):
        ptsUp = wing.airfoils[idx].ptsUp
        ptsLo = wing.airfoils[idx].ptsLo
        apex = wing.secApex[idx]
        chord = wing.chords[idx]
        angle = wing.secAngles[idx]
        return ptsUp, ptsLo, apex, chord, angle

    def write_curve(fid1, nConnect, pts, apex, chord, angle, rotAxis=0.25, rev=False):
        if rev:
            pts = np.flipud(pts)
        apex = np.array([apex[0],apex[2],apex[1]])
        pts = rotate_2d(pts, [rotAxis,0], angle) *chord
        pts = np.hstack([pts,np.zeros([len(pts),1])])
        pts += apex
        for pt in pts:
            fid1.write('  $_TMP(PW_%d) addPoint {%.10f %.10f %.10f}\n'%(nConnect,pt[0], pt[1], pt[2]))
        

    fid = open(glfPath,'wt')
    fid.write('# Pointwise automation by Maxim Tyan 2014\n')
    fid.write('package require PWI_Glyph 2.17.2\n')
    
    fid.write('pw::Application setUndoMaximumLevels 5\n')
    fid.write('pw::Application reset\n')
    fid.write('pw::Application markUndoLevel {Journal Reset}\n')
    
    fid.write('pw::Application clearModified\n')
    
    fid.write('set _TMP(mode_1) [pw::Application begin Create]\n')
    fid.write('  set _TMP(PW_1) [pw::SegmentSpline create]\n')
    # 1st segment
    ptsUp, ptsLo, apex, chord, angle = get_section_data(ac.wing, 0)
    write_curve(fid, 1, ptsUp, apex, chord, angle,rev=True)
    fid.write('  $_TMP(PW_1) setSlope Akima\n')
    fid.write('  set _TMP(curve_1) [pw::Curve create]\n')
    fid.write('  $_TMP(curve_1) addSegment $_TMP(PW_1)\n')
    fid.write('  unset _TMP(PW_1)\n')
    fid.write('$_TMP(mode_1) end\n')
    fid.write('unset _TMP(mode_1)\n')
    fid.write('pw::Application markUndoLevel {Create Curve}\n')
    
    fid.write('unset _TMP(curve_1)\n')
    fid.write('set _TMP(mode_2) [pw::Application begin Create]\n')
    fid.write('  set _TMP(PW_2) [pw::SegmentSpline create]\n')
    fid.write('  set _DB(1) [pw::DatabaseEntity getByName "curve-1"]\n')
    fid.write('  $_TMP(PW_2) addPoint [list 1 0 $_DB(1)]\n')
    write_curve(fid, 2, ptsLo[1:-1], apex, chord, angle)
    fid.write(' $_TMP(PW_2) addPoint [list 0 0 $_DB(1)]\n')
    fid.write('  $_TMP(PW_2) setSlope Akima\n')
    fid.write('  set _TMP(curve_1) [pw::Curve create]\n')
    fid.write('  $_TMP(curve_1) addSegment $_TMP(PW_2)\n')
    fid.write('  unset _TMP(PW_2)\n')
    fid.write('$_TMP(mode_2) end\n')
    fid.write('unset _TMP(mode_2)\n')
    fid.write('pw::Application markUndoLevel {Create Curve}\n')
    # 2nd segment
    fid.write('unset _TMP(curve_1)\n')
    fid.write('set _TMP(mode_3) [pw::Application begin Create]\n')
    fid.write('  set _TMP(PW_3) [pw::SegmentSpline create]\n')
    ptsUp, ptsLo, apex, chord, angle = get_section_data(ac.wing, 1)
    write_curve(fid, 3, ptsUp, apex, chord, angle)
    fid.write('  $_TMP(PW_3) setSlope Akima\n')
    fid.write('  set _TMP(curve_1) [pw::Curve create]\n')
    fid.write('  $_TMP(curve_1) addSegment $_TMP(PW_3)\n')
    fid.write('  unset _TMP(PW_3)\n')
    fid.write('$_TMP(mode_3) end\n')
    fid.write('unset _TMP(mode_3)\n')
    fid.write('pw::Application markUndoLevel {Create Curve}\n')
    
    fid.write('unset _TMP(curve_1)\n')
    fid.write('set _TMP(mode_4) [pw::Application begin Create]\n')
    fid.write('  set _TMP(PW_4) [pw::SegmentSpline create]\n')
    fid.write('  set _DB(2) [pw::DatabaseEntity getByName "curve-3"]\n')
    fid.write('  $_TMP(PW_4) addPoint [list 0 0 $_DB(2)]\n')
    write_curve(fid, 4, ptsLo[1:-1], apex, chord, angle)
    
    fid.write('  $_TMP(PW_4) addPoint [list 1 0 $_DB(2)]\n')
    fid.write('  $_TMP(PW_4) setSlope Akima\n')
    fid.write('  set _TMP(curve_1) [pw::Curve create]\n')
    fid.write('  $_TMP(curve_1) addSegment $_TMP(PW_4)\n')
    fid.write('  unset _TMP(PW_4)\n')
    fid.write('$_TMP(mode_4) end\n')
    fid.write('unset _TMP(mode_4)\n')
    fid.write('pw::Application markUndoLevel {Create Curve}\n')
    
    fid.write('unset _TMP(curve_1)\n')
    fid.write(' set _TMP(mode_5) [pw::Application begin Create]\n')
    fid.write('  set _TMP(PW_5) [pw::SegmentSpline create]\n')
    # tip
    ptsUp, ptsLo, apex, chord, angle = get_section_data(ac.wing, 2)
    write_curve(fid, 5, ptsUp, apex, chord, angle)

    fid.write('  $_TMP(PW_5) setSlope Akima\n')
    fid.write('  set _TMP(curve_1) [pw::Curve create]\n')
    fid.write('  $_TMP(curve_1) addSegment $_TMP(PW_5)\n')
    fid.write('  unset _TMP(PW_5)\n')
    fid.write('$_TMP(mode_5) end\n')
    fid.write('unset _TMP(mode_5)\n')
    fid.write('pw::Application markUndoLevel {Create Curve}\n')
    
    fid.write('unset _TMP(curve_1)\n')
    fid.write('set _TMP(mode_6) [pw::Application begin Create]\n')
    fid.write('  set _TMP(PW_6) [pw::SegmentSpline create]\n')
    fid.write('  set _DB(3) [pw::DatabaseEntity getByName "curve-5"]\n')
    fid.write('  $_TMP(PW_6) addPoint [list 1 0 $_DB(3)]\n')
    write_curve(fid, 6, ptsLo[1:-1], apex, chord, angle,rev=True)
    fid.write('$_TMP(PW_6) addPoint [list 0 0 $_DB(3)]\n')
    fid.write('$_TMP(PW_6) setSlope Akima\n')
    fid.write('set _TMP(curve_1) [pw::Curve create]\n')
    fid.write('$_TMP(curve_1) addSegment $_TMP(PW_6)\n')
    fid.write('unset _TMP(PW_6)\n')
    fid.write('$_TMP(mode_6) end\n')
    fid.write('unset _TMP(mode_6)\n')
    fid.write('pw::Application markUndoLevel {Create Curve}\n')
    
    fid2 = open('./templates/fw_db_automation_surface.txt','rt')
    for line in fid2:
        fid.write(line)
    fid2.close()
    fid.close()



def create_fw_cmesh(ac=None,igsPath=None,casPathSym=None, casPathNonsym=None, yplus=0.5,
                    glfPath=None):
    if ac==None:
        ac = aircraft_FW.load('Baseline1')
    if igsPath==None:
        igsPath = r'E:/1. Current work/2014 - UAV performance code/CFD automation/fw3.igs'
    if casPathSym==None:
        casPathSym = 'fw_case_sym.cas'
    if casPathNonsym==None:
        casPathNonsym = 'fw_case_nonsym.cas'
    if glfPath==None:
        glfPath='dbg_cmesh.glf'
    
    create_wing_db(ac,glfPath)

    fid = open('./templates/fw_automation.glf')
    lines = fid.readlines()
    fid.close()
    #print ac.get_cg()

    nWallChord = 65
    nWallSeg1  = 30
    nWallSeg2  = 31
    nWallTip   = 8
    nFFx = 46
    nFFy = 75
    nFFz = 46
    
    #yplus = 0.5
    dsLE = 1e-3

    dsFFy = 1e-3
    dsTopAft = 1e-2
    dsTopFront = 1e-2
    dsZ = 1e-3
    
    # igs input file
    
    #lines[14] = '  $_TMP(mode_1) initialize -type Automatic {%s}\n'%igsPath
    # wing tip
    tc = ac.wing.airfoils[-1].thickness
    tc *= ac.wing.chords[-1] / 2.
    sweep = ac.wing.segSweepLErad[-1] + np.radians(5)
    pt1offset = np.array([tc*np.tan(sweep), 0, tc])
    pt2offset = np.array([0, 0, tc])
    
    create_wing_db(ac,glfPath)
    
    dl = 34
    lines[40-dl] = '  $_TMP(PW_2) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 0]] {%.6f %.6f %.6f}]\n'%(pt1offset[0], pt1offset[1],pt1offset[2])
    lines[53-dl] = '  $_TMP(PW_3) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 1]] {%.6f %.6f %.6f}]\n'%(pt2offset[0], pt2offset[1],pt2offset[2])
    
    # ff aft
    span = ac.wing.span + tc
    lines[186-dl] = '  $_TMP(PW_13) addPoint [pwu::Vector3 add [$_CN(11) getPosition -arc 1] {0 0 %.6f}]\n'%span
    lines[198-dl] = '  $_TMP(PW_14) addPoint [pwu::Vector3 add [$_CN(13) getPosition -arc 1] {0 0 %.6f}]\n'%span
    lines[211-dl] = '  $_TMP(PW_15) addPoint [pwu::Vector3 add [$_CN(12) getPosition -arc 1] {0 0 %.6f}]\n'%span
    
    # number of points
    lines[380-dl] = '$_TMP(PW_28) do setDimension %d\n'%nWallChord
    lines[387-dl] = '$_TMP(PW_29) do setDimension %d\n'%nWallSeg1
    lines[394-dl] = '$_TMP(PW_30) do setDimension %d\n'%nWallSeg2
    lines[401-dl] = '$_TMP(PW_31) do setDimension %d\n'%nWallTip
    lines[408-dl] = '$_TMP(PW_32) do setDimension %d\n'%(nWallSeg1+nWallSeg2+nWallTip-2)
    
    lines[416-dl] = '$_TMP(PW_33) do setDimension %d\n'%nFFx
    lines[424-dl] = '$_TMP(PW_34) do setDimension %d\n'%nFFy
    lines[432-dl] = '$_TMP(PW_35) do setDimension %d\n'%nFFz
    lines[439-dl] = '$_TMP(PW_36) do setDimension %d\n'%(2*nWallChord-1)
    lines[446-dl] = '$_TMP(PW_37) do setDimension %d\n'%nFFz
    # wall spacing BL
    ds1 = ac.designGoals.fc.get_wall_spacing(yplus)
    cr = ac.wing.chords[0]
    ct = ac.wing.chords[-1]
    
    #FIXME: high taper ratio caused negative jacobean, so it is decided
    # temporarily assume yplus at MAC for all wing
    cr = ac.wing.MAC
    ct = ac.wing.MAC
    lines[453-dl] = '  $_TMP(PW_38) setBeginSpacing %.6e\n'%(ds1*cr)
    lines[456-dl] = '  $_TMP(PW_39) setBeginSpacing %.6e\n'%(ds1*cr)
    lines[464-dl] = '  $_TMP(PW_40) setBeginSpacing %.6e\n'%(ds1*ct)
    lines[467-dl] = '  $_TMP(PW_41) setBeginSpacing %.6e\n'%(ds1*ct)

    # leading/trailing edge spacing
    c2 = ac.wing.chords[1]
    dsLEr = dsLE*cr
    dsLE2 = dsLE*c2
    dsLEt = dsLE*ct
    
    lines[551-dl] = '  $_TMP(PW_59) setBeginSpacing %.10f\n'%(dsLEr)
    lines[554-dl] = '  $_TMP(PW_60) setEndSpacing %.10f\n'%(dsLEr)
    lines[557-dl] = '  $_TMP(PW_61) setBeginSpacing %.10f\n'%(dsLEr)
    lines[560-dl] = '  $_TMP(PW_62) setEndSpacing %.10f\n'%(dsLEr)
    
    lines[568-dl] = '  $_TMP(PW_63) setBeginSpacing %.10f\n'%(dsLE2)
    lines[571-dl] = '  $_TMP(PW_64) setEndSpacing %.10f\n'%(dsLE2)
    lines[574-dl] = '  $_TMP(PW_65) setBeginSpacing %.10f\n'%(dsLE2)
    lines[577-dl] = '  $_TMP(PW_66) setEndSpacing %.10f\n'%(dsLE2)
    
    lines[585-dl] = '  $_TMP(PW_67) setBeginSpacing %.10f\n'%(dsLEt)
    lines[588-dl] = '  $_TMP(PW_68) setEndSpacing %.10f\n'%(dsLEt)
    lines[591-dl] = '  $_TMP(PW_69) setBeginSpacing %.10f\n'%(dsLEt)
    lines[594-dl] = '  $_TMP(PW_70) setEndSpacing %.10f\n'%(dsLEt)
    lines[597-dl] = '  $_TMP(PW_71) setBeginSpacing %.10f\n'%(dsLEt)
    lines[600-dl] = '  $_TMP(PW_72) setEndSpacing %.10f\n'%(dsLEt)
    lines[608-dl] = '  $_TMP(PW_73) setEndSpacing %.10f\n'%(dsLEt)

    # case file export
    lines[1089] = '  $_TMP(mode_1) initialize -type CAE {%s}\n'%casPathNonsym
    lines[1113] = '  $_TMP(mode_4) initialize -type CAE {%s}\n'%casPathSym

    #glfPath = os.path.abspath('temp/new.glf')
    fid = open(glfPath,'at')
    for line in lines:
        fid.write(line)
    fid.close()
    #print glfPath
    os.system('\"%s\"'%glfPath)
    #os.remove(glfPath)


if __name__=="__main__":
    create_fw_cmesh()