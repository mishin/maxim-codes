# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 15:42:05 2014

@author: Maxim
"""
import convert
from paths import MyPaths
import os
from misc_tools import read_tabulated_data_without_header

pth = MyPaths()

def get_parasite_drag_fw(ac,altitude=0.0):
    
    path1 = pth.aeroCD0in
    path2 = pth.aeroCD0inWave
    
    W = ac.designGoals.grossMass*9.81
    rho = ac.designGoals.fc.atm.density
    V = ac.designGoals.cruiseSpeed
    Sref = ac.wing.area
    CL=W/(0.5*rho*V*V*Sref)

    # --- friction + form drag input ---
    fid1 = open(path1,'wt')
    #h = ac.designGoals.cruiseAltitude
    #h = convert.m_to_ft(h)
    h = float(altitude)
    Sref = convert.sqm_to_sqft(Sref)
    Abase = 0.095
    Dexit = 0.0
    Ewdd = 1.8
    prleak = 10.0
    fid1.write('%d %.0f %.0f %.4f %.1f %.1f %.1f\n'%(1,h,Sref,Abase, Dexit, Ewdd, prleak))

    ale = ac.wing.equivSweepLEdeg
    aqc = ac.wing.equivSweepC4deg
    ahc = ac.wing.equivSweepC2deg
    a4 = ac.wing.span
    a4 = convert.m_to_ft(a4)
    # ---
    t1 = ac.wing.airfoils[0].thickness*ac.wing.chords[0]
    t2 = ac.wing.airfoils[1].thickness*ac.wing.chords[1]
    Amax = (t1+t2)*ac.wing.segSpans[0]
    # ---
    ttc = ac.wing.equivThickness
    fc = ac.wing.equivCamber
    cld = CL
    xa = 1.02
    fid1.write('%.4f %.4f %.4f %.4f %.4f %.4f %.6f %.4f %.2f\n'%(ale,aqc,ahc,a4,Amax,ttc,fc,cld,xa))

    wSwet = ac.wing.wettedArea
    wSwet = convert.sqm_to_sqft(wSwet)
    wRefL = convert.m_to_ft(ac.wing.MAC)

    fid1.write('Wing %.4f %.4f 1 '%(wSwet,wRefL))
    wtc = ac.wing.equivThickness
    wc4 = ac.wing.equivSweepC4deg
    wxtc = ac.wing.equivThicknessLoc
    q = 1.0 #???
    fid1.write('%.4f %.4f %.4f %.4f'%(wtc,wc4,wxtc,q))
    fid1.close()

    # --- wave drag input ---
    fid2 = open(path2,'wt')
    fid2.write('1 0 0 0 0 0 0 0\n')
    wrle = ac.wing.equivLEradius
    cam = ac.wing.equivCamber
    fid2.write('%.4f %.4f %.4f %.4f %d\n'%(wtc,wxtc,wrle,cam, 1))
    fid2.write('0.0394 0.5 0. 0. 0\n')
    fid2.write('0. 0. 0. 0. 1\n')
    fid2.write('0.04 0.5 0. 0. 0\n')
    wte = ac.wing.equivSweepTEdeg
    lref = convert.m_to_ft(ac.wing.MAC)
    AR = ac.wing.aspectRatio
    fid2.write('%.4f %.4f %.4f %.4f %.4f\n'%(ale,wte,lref,AR,Sref))
    fid2.write('44.56 0. 0.29 2.24\n')
    fid2.write('0. 0. 0. 0.\n')
    fid2.write('44. 3. 0.277 2.48\n')
    fid2.write('0. 0. 0. 0. 0.\n')
    fid2.write('2\n')
    fid2.write('21.08 20.92 19.364 4.97 2.9\n')
    fid2.write('0.35 0.45 0.47 0.9339\n')
    fid2.close()
    os.system(pth.aeroCD0)

    data = read_tabulated_data_without_header(pth.aeroCD0out,3)
    pth.clean_drag_files()
    CD = data[:-1,-1]
    M = data[:-1,0]
    Mdd = data[-1,0]
    CDdd = data[-1,-1]
    return M, CD, Mdd, CDdd


def run_test1():
    import aircraft_FW
    ac = aircraft_FW.load('X-1')
    print get_parasite_drag_fw(ac)


if __name__=="__main__":
    run_test1()