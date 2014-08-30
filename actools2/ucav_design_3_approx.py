# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 11:54:45 2014

@author: Maxim
"""

from misc_tools import read_tabulated_data_without_header, Normalization, RbfMod
from scipy.optimize import fmin_slsqp, minimize
import numpy as np
from scipy.interpolate import Rbf

from ucav_design_1 import DesignFormulation

def get_rbf(x,f):
    f = np.vstack([x,f])
    f = np.transpose(f)
    return ModRbf(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8])


def run_optimization_1():
    pathIn = 'design_out4.txt'
    pathInSamples = 'DOE_LHS250_FFD2_3_3.txt'
    learnData = read_tabulated_data_without_header(pathIn)
    xNorm = read_tabulated_data_without_header(pathInSamples)
    #xNorm = np.transpose(xNorm)
    learnData = np.transpose(learnData)
    
    WemptyMax      = 3400.0
    CnbMin         = 0.0002
    ClbMax         = -0.05
    SMmin          = -0.05
    SMmax          = 0.10
    combatRadiusMin= 1000.0
    RCmin          = 125.0
    VmaxMin        = 0.9 # Mach
    
    LD      = RbfMod(xNorm, learnData[0])
    We      = RbfMod(xNorm, learnData[1])
    _Cnb     = RbfMod(xNorm, (learnData[2])*1e3)
    Clb     = RbfMod(xNorm, learnData[3])
    _SM      = RbfMod(xNorm, learnData[4]*1e3)
    Rcombat = RbfMod(xNorm, learnData[5])
    RC      = RbfMod(xNorm, learnData[6])
    Vmax    = RbfMod(xNorm, learnData[7])
    
    Cnb = lambda x: _Cnb(x)/1e3
    SM = lambda x: _SM(x)/1e3
    bnds = np.ones([9,2])
    bnds[:,0] = -bnds[:,0]
    bnds = tuple(bnds)

    f  = lambda x: -LD(x)
    g1 = lambda x: WemptyMax - We(x)
    g2 = lambda x: (Cnb(x) - CnbMin)*1e3
    g3 = lambda x: ClbMax - Clb(x)
    g4 = lambda x: SMmax - SM(x)
    g5 = lambda x: SM(x) - SMmin
    g6 = lambda x: Rcombat(x)/1e3 - combatRadiusMin
    g7 = lambda x: RC(x) - RCmin
    g8 = lambda x: Vmax(x) - VmaxMin

    cnstr = ({'type':'ineq', 'fun':g1},
             {'type':'ineq', 'fun':g2},
             {'type':'ineq', 'fun':g3},
             {'type':'ineq', 'fun':g4},
             {'type':'ineq', 'fun':g5},
             {'type':'ineq', 'fun':g6},
             {'type':'ineq', 'fun':g7},
             {'type':'ineq', 'fun':g8})

    x0 = np.zeros(9)
    rslt = minimize(f, x0, method='SLSQP',constraints=cnstr, bounds=bnds)
    xopt = rslt.x
    print rslt
    print g1(xopt), We(xopt)
    print g2(xopt), Cnb(xopt)
    print g3(xopt), Clb(xopt)
    print g4(xopt), SM(xopt)
    print g5(xopt), SM(xopt)
    print g6(xopt), Rcombat(xopt)
    print g7(xopt), RC(xopt)
    print g8(xopt), Vmax(xopt)

    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    normalize = Normalization(ac.lb,ac.ub)
    xoptDenorm = normalize.denormalize(xopt)
    print xoptDenorm
    ac.set_x(xoptDenorm)
    print ac.analysisData, '\n---'
    print (ac.analysisData[0] - LD(xopt)) / ac.analysisData[0]
    print (ac.analysisData[1] - We(xopt)) / ac.analysisData[1]
    print (ac.analysisData[2] - Cnb(xopt)) / ac.analysisData[2]
    print (ac.analysisData[3] - Clb(xopt)) / ac.analysisData[3]
    print (ac.analysisData[4] - SM(xopt)) / ac.analysisData[4]
    print (ac.analysisData[5] - Rcombat(xopt)) / ac.analysisData[5]
    print (ac.analysisData[6] - RC(xopt)) / ac.analysisData[6]
    print (ac.analysisData[7] - Vmax(xopt)) / ac.analysisData[7]

    ac.display()


if __name__=="__main__":
    run_optimization_1()