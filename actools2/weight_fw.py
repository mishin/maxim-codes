# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 21:39:57 2014

@author: Maxim
"""
import numpy as np
import constants
from weight_tools import MassList, MassComponent, AircraftMass, Fuel
import convert


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
        mEnM, cgEnM = self._mass_engine_mount()
        mEnS, cgEnS = self._mass_engine_section()
        mEnC, cgEnC = self._mass_engine_cooling()
        mAir, cgAir = self._mass_airintake()
        mOil, cgOil = self._mass_oil_cooling()
        self.total.airframe.add_item('wing',mWing,cgWing)
        self.total.airframe.add_item('fuselage',mFus,cgFus)
        self.total.airframe.add_item('Main landing gear',mMlg,cgMlg)
        self.total.airframe.add_item('Nose landing gear',mNlg,cgNlg)
        self.total.airframe.add_item('Engine mount',mEnM,cgEnM)
        self.total.airframe.add_item('Engine section',mEnS,cgEnS)
        self.total.airframe.add_item('Engine cooling',mEnC,cgEnC)
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
        self.Wdg = convert.kg_to_lb(self.ac.designGoals.grossMass)
        self.Nl = self.ac.designGoals.loadFactorLanding
        self.Wt = self.Wdg
        self.Wl = self.Wdg
        self.Ne = self.ac.propulsion.numberOfEngines
        self.T  = self.ac.propulsion.totalThrust
        
    def _mass_wing(self,wing):
        # by Raymer - fighter weights
        Sw = convert.sqm_to_sqft(wing.area)
        A  = wing.aspectRatio
        tcRoot = wing.airfoils[0].thickness
        lmbda = wing.sweepElasticRad
        Scsw  = 0.1*Sw #FIXME: calculate control surface area
        mWing = 0.0103*self.Kdw*self.Kvs* (self.Wdg*self.Nz)**0.5* Sw**0.622* A**0.785 *tcRoot**(-0.4)* (1+lmbda)**0.5*np.cos(lmbda)**(-1.0)* Scsw**0.04
        return convert.lb_to_kg(mWing), np.zeros(3)
    
    def _mass_fuselage(self):
        L = convert.m_to_ft(self.ac.wing.chords[0])
        D = convert.m_to_ft(self.ac.wing.airfoils[0].thickness)
        W = convert.m_to_ft(self.ac.fusWidth)
        mFus = 0.499*self.Kdwf*self.Wdg**0.35 *self.Nz**0.25 *L**0.5 *D**0.849 *W**0.685
        return convert.lb_to_kg(mFus), np.zeros(3)
    
    def _mass_mlg(self):
        Lm = convert.m_to_ft(self.ac.landingGear.strutLength[1])
        mMlg = self.Kcb*self.Ktpg*(self.Wt*self.Nl)**0.25 * Lm**0.973
        return convert.lb_to_kg(mMlg), np.zeros(3)

    def _mass_nlg(self):
        Ln = convert.m_to_ft(self.ac.landingGear.strutLength[0])
        Nnw = 1.0 # number of nose wheels
        mNlg = (self.Wl*self.Nl)**0.29 *Ln**0.5 *Nnw**0.525
        return convert.lb_to_kg(mNlg), np.zeros(3)

    def _mass_engine_mount(self):
        mEngineMount = 0.013*self.Ne**0.795*self.T**0.579 *self.Nz
        return convert.kg_to_lb(mEngineMount), np.zeros(3)
        
    def _mass_engine_section(self):
        We = convert.kg_to_lb(self.ac.propulsion.engineMass)
        mEngineSection = 0.01*We**0.717*self.Ne*self.Nz
        return convert.lb_to_kg(mEngineSection), np.zeros(3)
        
    def _mass_airintake(self):
        #mAir = 13.29*Kvg * Ld**0.643* Kd**0.182 *Nen**1.498 * (Ls/Ld)**(-.373)*De
        return 0.0, np.zeros(3)
        
    def _mass_engine_cooling(self):
        return 0.0, np.zeros(3)
    def _mass_oil_cooling(self):
        return 0.0, np.zeros(3)
    