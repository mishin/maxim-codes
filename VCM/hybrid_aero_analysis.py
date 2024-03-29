# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 11:07:19 2014

@author: Maxim
"""

import airfoil
import CFDpolar as cfd
import numpy as np
from fluent_solver import FluentOutput
from scipy.interpolate import interp1d
import os
from miscTools import denormalize
from scipy.optimize import fminbound
from misc_tools import read_tabulated_data_without_header

def write_polar(path,polar,polarType):
    fid = open(path,'wt')
    fid.write('%s\n'%polarType)
    fid.write('alpha\tcl\tcd\tcm\n')
    for alpha, cl, cd, cm in zip(polar.alpha, polar.cl, polar.cd, polar.cm):
        fid.write('%.2f\t%.6f\t%.6f\t%.6f\n'%(alpha,cl,cd,cm))
    fid.close()

def aero_analysis(af, fc, filename=None,alpha=None):
    if alpha==None:
        alpha = np.arange(0,17,2.0)
    result1 = af.get_J_polar(fc.Mach, fc.Re, [alpha[0],alpha[-1],0.5])
    solver = cfd.CFDsolver(af,fc,1.0,mesh='O')
    solver.fluent.residuals['energy']=1e-6
    solver.fluent.relaxationFactor['xvelocity'] = 1e-3
    solver.mesh._airfoilPts = 75
    solver.mesh._interiorPts = 75
    solver.mesh._dsTE = 5e-5
    solver.mesh._dsLE = 2e-3
    solver.mesh._growthRate = 1.2
    solver.create_mesh()
    
    result2 = solver.run_for_multiple_aoa(alpha,turbulenceModel='ke-realizable')
    result3 = combine_results(result1, result2)
    result3.thickness = af.thickness
    if not filename==None:
        write_polar(filename, result3,'combined')
    return result3


def combine_results(resultXfoil, resultCfd):
    alphaMax = min([max(resultXfoil.alpha), max(resultCfd.alpha)])
    alphaMin = max([min(resultXfoil.alpha), min(resultCfd.alpha)])
    alphaNew = np.arange(alphaMin,alphaMax+1.,1.)
    clAlpha = interp1d(resultCfd.alpha,resultCfd.cl,'cubic')
    cdAlpha = interp1d(resultXfoil.alpha, resultXfoil.cd,'cubic')
    cmAlpha = interp1d(resultXfoil.alpha, resultXfoil.cm,'cubic')
    resultNew = FluentOutput()
    resultNew.alpha = alphaNew
    resultNew.cl = clAlpha(alphaNew)
    resultNew.cd = cdAlpha(alphaNew)
    resultNew.cm = cmAlpha(alphaNew)
    f1 = lambda x: -clAlpha(x)
    f2 = lambda x: -clAlpha(x)/cdAlpha(x)
    f3 = lambda x: -clAlpha(x)**1.5/cdAlpha(x)
    a1 = fminbound(f1,alphaMin, alphaMax, full_output=1)
    a2 = fminbound(f2,alphaMin, alphaMax, full_output=1)
    a3 = fminbound(f3,alphaMin, alphaMax, full_output=1)
    resultNew.alphaClmax = a1[0]
    resultNew.alphaLDmax = a2[0]
    resultNew.alphaLD32max = a3[0]
    resultNew.Clmax = -a1[1]
    resultNew.LDmax = -a2[1]
    resultNew.LD32max = -a3[1]
    resultNew.cdAtLDmax = cdAlpha(a2[0])
    resultNew.cdAtLD32max = cdAlpha(a3[0])
    return resultNew

def run_test1():
    A = np.array([0.15225392, 0.17285363, 0.15722004, 0.13938069,-0.08884681, -0.05170833, -0.01974148, 0.02985025])
#    A[3] += -0.04
#    A[7] += 0.04
    af = airfoil.cst(A[0:4],A[4:8])
    af.plot()
    fc = cfd.FlightConditions(22.0,0.0,0,0.24)
    result = aero_analysis(af,fc,'testPolar.txt')

def calculate_bounds():
    X = np.array([0.14391813, 0.18778261, 0.14634264, 0.15348147, 0.15107265, -0.09014438, -0.05862712, -0.03488944, -0.01428362, 0.03831908])
    lb = X - np.array([0.01, 0.05, 0.05, 0.05, 0.00, 0.05, 0.05, 0.05, 0.05, 0.05])
    ub = X + np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.01, 0.05, 0.05, 0.05, 0.00])
    af1 = airfoil.cst(lb[0:5],ub[5:10])
    af1.plot()
    af2 = airfoil.cst(ub[0:5],lb[5:10])
    af2.plot()
    
    fc = cfd.FlightConditions(22.0,0.0,0,0.24)
    result = aero_analysis(af1,fc,'thin.txt')
    result = aero_analysis(af2,fc,'thick.txt')


def read_doe(path):
    fid = open(path,'rt')
    lines = fid.readlines()
    output = np.zeros([100,10])
    fid.close()
    i = 0
    for line in lines:
        if not line.strip()=='':
            seg = line.split()
            for j,val in enumerate(seg):
                output[i,j] = float(val)
            i += 1
    return output[:i]

def write_output(x,result, path,clCruise):
    fid = open(path,'a')
    for val in x:
        fid.write('%.6f\t'%val)
    fid.write('%.6f\t'%result.LDmax)
    fid.write('%.6f\t'%result.alphaLDmax)
    fid.write('%.6f\t'%result.Clmax)
    fid.write('%.6f\t'%result.alphaClmax)
    fid.write('%.6f\t'%result.cdAtLDmax)
    fid.write('%.6f\t'%result.LD32max)
    fid.write('%.6f\t'%result.thickness)
    fid.write('%.6f\n'%result.get_cd_at_cl(clCruise))
    fid.close()

def get_initial_data():
    X = np.array([0.14391813, 0.18778261, 0.14634264, 0.15348147, 0.15107265, -0.09014438, -0.05862712, -0.03488944, -0.01428362, 0.03831908])
    lb = X - np.array([0.02, 0.02, 0.02, 0.02, 0.00, 0.02, 0.02, 0.02, 0.02, 0.02])
    ub = X + np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.00])
    fc = cfd.FlightConditions(22.0,0.0,0,0.24)
    return X,lb,ub,fc

def run_full_table_analysis(wdir):
    X,lb,ub,fc = get_initial_data()
    
    pathDOE = os.getcwd() + '\\' + wdir + '\\1. input.txt'
    pathOutput = os.getcwd() + '\\' + wdir + '\\2. outputAero.txt'
    clcruise = 0.4
    
    fid = open(pathOutput,'wt')
    fid.write('Au0\tAu1\tAu2\tAu3\tAu4\tAl0\tAl1\tAl2\tAl3\tAl4\tLDmax\talphaLDmax\tClmax\talphaClmax\tCdAtLDmax\tLD32max\tthickness\tcdcruise\n')
    fid.close()
    #designs = read_doe(pathDOE)
    designs = read_tabulated_data_without_header(pathDOE)
    for i,design in enumerate(designs):
        design = denormalize(design,lb,ub)
        af = airfoil.cst(design[0:5],design[5:10])
        filename = wdir + 'polar_%d.txt'%(i+1)
        result = aero_analysis(af, fc, filename)
        write_output(design,result,pathOutput,clcruise)
        af.write_txt(wdir + 'coordinates_%d.txt'%(i+1))

def run_baseline():
    #X = np.array([0.14391813, 0.18778261, 0.14634264, 0.15348147, 0.15107265, -0.09014438, -0.05862712, -0.03488944, -0.01428362, 0.03831908])
    #X = np.array([0.1639, 0.2078, 0.1663, 0.1629, 0.1511, -0.1049, -0.0548, -0.0149, -0.0009, 0.0383])
    X = np.array([0.1439, 0.1878,0.1463,0.1535,0.1511,-0.0901,-0.0586,-0.0349,-0.0143,0.0383])
    fc = cfd.FlightConditions(22.0,0.0,0,0.24)
    fid = open('RENNbaseline.txt','wt')
    fid.write('Au0\tAu1\tAu2\tAu3\tAu4\tAl0\tAl1\tAl2\tAl3\tAl4\tLDmax\talphaLDmax\tClmax\talphaClmax\tCdAtLDmax\tthickness\n')
    fid.close()
    #designs = read_doe(pathDOE)
    #for i,design in enumerate(designs):
    af = airfoil.cst(X[0:5],X[5:10])
    result = aero_analysis(af, fc, 'tmp1z.txt')
    write_output(X,result,'RENNbaseline.txt')

def run_multiple_tables():
    wdir = ['RENNaero\\75samples']
    for directory in wdir:
        run_full_table_analysis(directory)

if __name__=="__main__":
    run_baseline()