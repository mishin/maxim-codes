# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 17:07:48 2014

@author: Maxim
"""


import numpy as np
import db_tools
import engine
import propeller
import airfoil
from paths import MyPaths
import misc_tools

import drag
from mass import massList
from matplotlib.mlab import find
import AircraftAnalysis
import geometry as geom
from scipy import interp

class Aircraft(object):
    def __init__(self,dbPathAircraft):
        self.db = dbPathAircraft

    def load(self,name):
        pass
    def save(self,name):
        pass
    def display(self):
        pass
    
    def get_mass_empty(self,update=False):
        if update:
            self.update_mass()

    def get_mass_total(self,update=False):
        if update:
            self.update_mass()

    def get_cg(self,update=False):
        if update:
            self.update_mass()

    def get_inertia(self,update=False):
        if update:
            self.update_mass()

    def update_mass(self):
        pass

    def update_drag(self):
        pass

    def get_drag(self,update=False):
        if update:
            self.update_drag()
    
    def get_aero_trim(self,velocity=None,altitude=None,loadFactor=1.0,CmTrim=0.0,
                      flap=0.0):
        pass
    
    def get_aero(self,velocity=None,altitude=None,flap=0.0,elevator=0.0,
                 aileron=0.0,rudder=0.0):
        pass

class Wing(object):
    def __init__(self):
        self.dbPathAirfoil=paths.Database().airfoil
        self.segLengths   =None
        self.chords       =None
        self.airfoilNames =None
        self.segOffsets   =None
        self.segDihedrals =None
        self.twistAngles  =None
        self.incidence    =None
        self.aapex        =None
        self.material     =None
        self.symmetry     =None
    
    def get_wing_data(self):
        self._analyze_geometry()
        self._create_mesh()
    
    def _analyze_geometry():
        pass
    def _create_mesh():
        pass
    def set_fuel_tank():
        pass
    

class MainWing(Wing):
    def add_aileron(self):
        pass
    def add_flap(self):
        pass

class HorizontalStab(Wing):
    def add_elevator(self):
        pass

class VerticalStab(Wing):
    def add_rudder(self):
        pass

class ControlSurface(object):
    pass