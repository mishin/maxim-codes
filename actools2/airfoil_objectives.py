# -*- coding: utf-8 -*-
"""
Created on Sat Mar 01 00:12:45 2014

@author: Maxim
"""

from airfoil import *
from misc_tools import Normalization

class AirfoilObjective:
    def __init__(self,x0,lb,ub):
        self.x0 = x0
        self.af = None
        self.clCruise = [0.2,0.3]
        self.clmaxMin = 1.5
        self.tcMin = 0.135
        self.tcMax = 0.145
        self.norm = Normalization(lb,ub)
    
    def set_cst(self,x):
        pass
    
    def f(self,x):
        pass
    
    def g1high(self,x):
        pass
    def g1low(self,x):
        pass
    def g2(self,x):
        pass
    def g3(self,x):
        pass