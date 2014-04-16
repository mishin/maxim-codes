# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:07:43 2014

@author: Maxim
"""
import numpy as np


class TurbofanEngine(object):
    def __init__(self):
        self.thrustSL = 0.0
        self.sfc = 0.0
        self.CGx = None
        self.CGy = None
        self.CGz = None
        self.mass = 0.0