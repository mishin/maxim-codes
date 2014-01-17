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

def write_polar(path,polar,polarType):
    fid = open(path,'wt')
    fid.write('%s\n'%polarType)
    fid.write('alpha\tcl\tcd\tcm\n')
    for alpha, cl, cd, cm in zip(polar.alpha, polar.cl, polar.cd, polar.cm):
        fid.write('%.2f\t%.6f\t%.6f\t%.6f\n'%(alpha,cl,cd,cm))
    fid.close()

def aero_analysis(af, fc, filename):
    result1 = af.get_X_polar(fc.Mach, fc.Re, [0,16,2.0],nIter=100)
    solver = cfd.CFDsolver(af,fc,1.0,mesh='O')
    solver.fluent.residuals['energy']=1e-6
    solver.fluent.relaxationFactor['xvelocity'] = 1e-3
    solver.mesh._airfoilPts = 75
    solver.mesh._interiorPts = 75
    solver.mesh._dsTE = 5e-5
    solver.mesh._dsLE = 2e-3
    solver.mesh._growthRate = 1.2
    solver.create_mesh()
    alpha = np.arange(0,17,2.0)
    result2 = solver.run_for_multiple_aoa(alpha,turbulenceModel='ke-realizable')
#    write_polar('xfoil_%s'%filename, result1, 'xfoil')
#    write_polar('cfd_%s'%filename, result2, 'cfd')
    result3 = combine_results(result1, result2)
    write_polar('combined_%s'%filename, result3,'combined')


def combine_results(resultXfoil, resultCfd):
    alphaMax = min([max(resultXfoil.alpha), max(resultCfd.alpha)])
    alphaMin = max([min(resultXfoil.alpha), min(resultCfd.alpha)])
    alphaNew = np.arange(alphaMin,alphaMax+1.,2.)
    clAlpha = interp1d(resultCfd.alpha,resultCfd.cl,'cubic')
    cdAlpha = interp1d(resultXfoil.alpha, resultXfoil.cd,'cubic')
    cmAlpha = interp1d(resultXfoil.alpha, resultXfoil.cm,'cubic')
    resultNew = FluentOutput()
    resultNew.alpha = alphaNew
    resultNew.cl = clAlpha(alphaNew)
    resultNew.cd = cdAlpha(alphaNew)
    resultNew.cm = cmAlpha(alphaNew)
    return resultNew

def run_test1():
    A = np.array([0.15225392, 0.17285363, 0.15722004, 0.13938069,-0.08884681, -0.05170833, -0.01974148, 0.02985025])
#    A[3] += -0.04
#    A[7] += 0.04
    af = airfoil.cst(A[0:4],A[4:8])
    af.plot()
    fc = cfd.FlightConditions(22.0,0.0,0,0.24)
    result = aero_analysis(af,fc,'testPolar.txt')


if __name__=="__main__":
    run_test1()