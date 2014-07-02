# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 15:42:05 2014

@author: Maxim
"""
import convert
from paths import MyPaths
from scipy.interpolate import Rbf
from flight_conditions import FlightConditions

pth = MyPaths()

def get_parasite_drag_fw():
    import aircraft_FW
    ac = aircraft_FW.load('X45C')
    
    W = ac.designGoals.grossMass*9.81
    rho = ac.designGoals.fc.atm.density
    V = ac.designGoals.cruiseSpeed
    Sref = ac.wing.area
    CL=W/(0.5*rho*V*V*Sref)
    
    path1 = 'input1.txt'
    path2 = 'input2.txt'
    fid1 = open(path1,'wt')
    fid2 = open(path2,'wt')

    h = ac.designGoals.cruiseAltitude
    h = convert.m_to_ft(h)
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
    wtc = ac.wing.equivThickness
    fid1.write('Wing %.4f %.4f %.0f %.4f %.4f %.4f 1.0\n'%(wSwet,wRefL,1,wtc,1))

def _get_parasite_drag_fw(aircraft,flightConditions):
    W = aircraft.designGoals.grossMass*9.81
    rho = flightConditions.atm.density
    V = flightConditions.velocity
    Sref = aircraft.wing.area
    CL=W/(0.5*rho*V**2*Sref)
    a1 = 1
    a2 = convert.m_to_ft(flightConditions.atm.altitude)
    a3 = convert.sqm_to_sqft(aircraft.wing.area)
    a4 = 0.095
    a5 = 0
    a6 = 1.8
    a7 = 10.0
    b1 = aircraft.wing.equivSweepLEdeg
    b2 = aircraft.wing.equivSweepC4deg
    b3 = aircraft.wing.equivSweepC2deg
    b4 = convert.m_to_ft(aircraft.wing.span)
    b5 = (aircraft.wing.airfoils[0].thickness+aircraft.wing.airfoils[1].thickness)*aircraft.wing.segSpans[0]
    b6 = aircraft.wing.airfoils[0].thickness
    b7 = aircraft.wing.equivCamber
    b8 = CL
    b9 = 1.02
    c1 = convert.sqm_to_sqft(aircraft.wing.wettedArea)
    c2 = convert.m_to_ft(aircraft.wing.span)
    c3 = 1
    c4 = aircraft.wing.equivThickness
    c5 = aircraft.wing.equivSweepC4rad
    c6 = aircraft.wing.equivThickness
    c7 = 1
    fid = open(pth.aeroCD0in,'wt')
    fid.write('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n'%(a1,a2,a3,a4,a5,a6,a7))
    fid.write('%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n'%(b1,b2,b3,b4,b5,b6,b7,b8,b9))
    fid.write('Wing\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f'%(c1,c2,c3,c4,c5,c6,c7))
    fid.close()
    # wave drag input
    w1 = aircraft.wing.equivThickness
    w2 = aircraft.wing.equivThicknessLoc
    w3 = aircraft.wing.equivLEradius
    w4 = aircraft.wing.equivCamber
    w5 = 1
    w11 = aircraft.wing.equivSweepLEdeg
    w12 = -aircraft.wing.equivSweepTEdeg
    w13 = aircraft.wing.taperRatio
    w14 = aircraft.wing.aspectRatio
    w15 = convert.sqm_to_sqft(aircraft.wing.area)
    fid2 = open(pth.aeroCD0inWave,'wt')
    fid2.write('1\t0\t0\t0\t0\t0\t0\t0\n')
    fid2.write('%.6f\t%.6f\t%.6f\t%.6f\t%d\n'%(w1,w2,w3,w4,w5))
    # fixed input for parts except wing. this input is not considered in calc
    fid2.write('0.155 0.5 1.58 0. 1.\n')
    fid2.write('0. 0. 0. 0. 0.\n')
    fid2.write('0.155 0.5 1.58 0. 1.\n')
    #
    fid2.write('%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n'%(w11,w12,w13,w14,w15))
    fid2.write('35. 13. 0.3 3.96\n')
    fid2.write('0. 0. 0. 0.\n')
    fid2.write('49. 13. 0.27 3.5\n')
    fid2.write('0. 0. 0. 0. 0.\n')
    fid2.write('2\n')
    fid2.write('6.2 1.43 1.17 3. 4.34\n')
    fid2.write('0. 0. 0. 0.')
    fid2.close()


def run_test1():
    import aircraft_FW
    from flight_conditions import FlightConditions
    ac = aircraft_FW.load('X47A')
    fc = FlightConditions(0.7,10000,0,ac.wing.MAC)
    get_parasite_drag_fw(ac,fc)

def run_test2():
    get_parasite_drag_fw()
    
if __name__=="__main__":
    run_test2()