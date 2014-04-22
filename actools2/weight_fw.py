# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 21:39:57 2014

@author: Maxim
"""
import numpy as np
import constants

from weight_tools import MassList, MassComponent, AircraftMass, Fuel

def get_flying_wing_mass(aircraft):
    mass = BlendedWingBodyMass(aircraft)
    mass.analyze()


class BlendedWingBodyMass(object):
    def __init__(self,aircraft):
        self.ac = aircraft
    def analyze(self):
        self._analyze_wing(self.ac.wing)
    def _analyze_wing(self,wing):
        #FIXME: assumptions that need to be fixed
        bcs = 0.5 # central (non bending) section
        etaCp = wing.macLocation[1]
        rhoMat = 2700.0
        # ----
        g = constants.GRAVITY_ACCEL
        nUlt = self.ac.designGoals.loadFactor
        b = wing.span
        Wg = self.ac.designGoals.grossWeight
        Rtip = 1.0 # if no winglets
        RL = 1.0 # for flying wing
        Rcs = (1.0-bcs/b)*(1+bcs/b*(2.45*np.cos(wing.sweepElacticRad)-1.0))
        bSt = wing.structuralSpan
        Wbox = 0.36*nUlt*Rtip*Rcs*RL*Wg*etaCp*bSt*rhoMat*g/sigmaBar*(1.05*Rcant/etaT+3.67)
        Wrib = rhoMat*g*krib*S*(tref+(tr+tt)/2.0)