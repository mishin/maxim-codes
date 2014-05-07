# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 13:06:22 2014

Main class for flying wing configuration aerodynamic analysis using AVL code.
This class can be used as parent for different configurations analysis.

1. code for flying wing analysis with elevator only

@author: Maxim
"""

from aero_avl_tools import AVLresults
from paths import MyPaths
import os
import numpy as np
import shlex
from StringIO import StringIO

pth = MyPaths()


class AVLsolver(object):
    def __init__(self,aircraft):
        self.ac = aircraft

    def create_input_file(self):
        """
        basic terms
        wing configuration
        control surface - ELEVATOR only for now
        weight and CG
        """
        pth.set_file_prefix('flyingwing')
        path = pth.get_tmp_file('avl')
        cg = self.ac.get_cg()
        CD0 = self.ac.get_drag()
        
        fid = open(path,'wt')
        # general input
        fid.write('%s\n'%self.ac.name)
        fid.write('0.0\n')
        fid.write('0  0  0.0\n')
        fid.write('%.4f  %.4f  %.4f\n'%(self.ac.wing.area,self.ac.wing.MAC,self.ac.wing.span))
        fid.write('%.4f  %.4f  %.4f\n'%(cg[0],cg[1],cg[2]))
        fid.write('%.8f\n'%CD0)
        # wing input
        fid.write('SURFACE\n')
        fid.write('WING\n')
        fid.write('%d  0  %d  0\n'%(self.ac.vlm.panelsChordwise, self.ac.vlm.panelsSpanwise))
        fid.write('YDUPLICATE\n0\n')
        for i in range(self.ac.wing.nSec):
            fid.write('#--- section %d ---\n'%(i+1))
            pathAf = pth.get_tmp_file('txt','af%d'%(i+1))
            self.ac.wing.airfoils[i].write_txt(pathAf)
            apex = self.ac.wing.secApex[i]
            chord = self.ac.wing.chords[i]
            angle = self.ac.wing.secAngles[i]
            fid.write('SECTION\n')
            fid.write('%.6f  %.6f  %.6f  %.6f  %.6f\n'%(apex[0],apex[1],apex[2],chord,angle))
            fid.write('AFIL\n%s\n'%pathAf)
            fid.write('CONTROL\n')
            fid.write('elevator  1.0  %.4f  0  0  0  +1\n'%(self.ac.wing.elevon.location[i]))
        fid.close()
    
    def run_trim(self):
        self.create_input_file()
        self.run_avl()
        self.collect_data()

    def set_flight_conditions(self,flightConditions):
        pass
    def run_single_point(self):
        # create input
        # run avl: load, oper single point, process data
        # clean files
        pass
    def run_alpha_sweep(self):
        pass
    def run_beta_sweep(self):
        pass


class Aerodynamics(AVLsolver):
    def __init__(self,aircraft):
        super(Aerodynamics,self).__init__(aircraft)

class RunCases:
    def __init__(self):
        pass

def run_test1():
    import aircraft_FW as aircraft
    ac = aircraft.load('sample_B45c')
    aero = Aerodynamics(ac)
    aero.create_input_file()

if __name__=="__main__":
    run_test1()