# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 13:06:22 2014

Main class for flying wing configuration aerodynamic analysis using AVL code.
This class can be used as parent for different configurations analysis.

@author: Maxim
"""

from aero_avl_tools import AVLresults
from paths import MyPaths
import os
import numpy as np
import shlex
from StringIO import StringIO

pth = MyPaths()

class RunCases:
    def __init__(self):
        pass

class AVLsolver(object):
    def __init__(self,aircraft):
        self.ac = aircraft
    def create_input_file(self,path):
        # flying wing configuration
        pass
    
    def run_trim(self):
        self.create_input_file(filepath)
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
