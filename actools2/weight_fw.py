# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 21:39:57 2014

@author: Maxim
"""
import numpy as np
import constants
import quantities as qu


from weight_tools import MassList, MassComponent, AircraftMass, Fuel

def get_flying_wing_mass(aircraft):
    mass = BlendedWingBodyMass(aircraft)
    mass.analyze()


class BlendedWingBodyMass(object):
    def __init__(self,aircraft):
        self.ac = aircraft
        self.total = AircraftMass(self.ac.name)
    def analyze(self):
        self._get_data()
        mWing, cgWing = self._mass_wing(self.ac.wing)
        mFus, cgFus = self._mass_fuselage()
        mMlg, cgMlg = self._mass_mlg()
        mNlg, cgNlg = self._mass_nlg()
        self.total.airframe.add_item('wing',mWing,cgWing)
        self.total.airframe.add_item('fuselage',mFus,cgFus)
        self.total.airframe.add_item('Main landing gear',mMlg,cgMlg)
        self.total.airframe.add_item('Nose landing gear',mNlg,cgNlg)
        self.total.display()
    
    def _get_coefficients(self):
        self.Kdw = 1.0 # 0.786 for delta wing, otherwise Kdw=1.0
        self.Kvs = 1.0 # 1.19 for variable sweep, otherwise 1.0
        self.Kdwf = 1.0 # 0.774 for delta wing, otherwise 1.0
        self.Kcb = 1.0
        self.Ktpg = 1.0

    def _get_data(self):
        self._get_coefficients()
        self.Nz = 1.5 * self.ac.designGoals.loadFactor # ultimate load factor
        self.Wdg = float((self.ac.designGoals.grossMass * qu.kg).rescale('pound'))
        
        
    def _mass_wing(self,wing):
        # by Raymer - fighter weights
        Sw = float((wing.area*qu.m**2).rescale('ft**2'))
        A  = wing.aspectRatio
        tcRoot = wing.airfoils[0].thickness
        lmbda = wing.sweepElasticRad
        Scsw  = 0.1*Sw #FIXME: calculate control surface area
        mWing = 0.0103*self.Kdw*self.Kvs* (self.Wdg*self.Nz)**0.5* Sw**0.622* A**0.785 *tcRoot**(-0.4)* (1+lmbda)**0.5*np.cos(lmbda)**(-1.0)* Scsw**0.04
        #Strapz = (wing.chords[0]+wing.chords[-1])*wing.span/2.0
        return float((mWing*qu.pound).rescale('kg')), np.zeros(3)
    
    def _mass_fuselage(self):
        L = float((self.ac.wing.chords[0]*qu.m).rescale('ft'))
        D = float((self.ac.wing.airfoils[0].thickness*qu.m).rescale('ft'))
        W = float((self.ac.fusWidth*qu.m).rescale('ft'))
        mFus = 0.499*self.Kdwf*self.Wdg**0.35 *self.Nz**0.25 *L**0.5 *D**0.849 *W**0.685
        return float((mFus*qu.pound).rescale('kg')), np.zeros(3)
    
    def _mass_mlg(self):
        Nl = self.ac.designGoals.loadFactorLanding
        Wt = self.Wdg
        Lm = float((self.ac.landingGear.strutLength[1]*qu.m).rescale('ft'))
        mMlg = self.Kcb*self.Ktpg*(Wt*Nl)**0.25 * Lm**0.973
        return float((mMlg*qu.pound).rescale('kg')), np.zeros(3)

    def _mass_nlg(self):
        return 0.0, np.zeros(3)
    def _mass_engine_mount(self):
        return 0.0, np.zeros(3)
    def _mass_engine_section(self):
        return 0.0, np.zeros(3)
    def _mass_airintake(self):
        return 0.0, np.zeros(3)
    def _mass_engine_cooling(self):
        return 0.0, np.zeros(3)
    def _mass_oil_cooling(self):
        return 0.0, np.zeros(3)
    