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
        self.fuelProp = constants.load('fuel_density')
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
        self.total.airframe.add_item('Oil cooling',mOil,cgOil)
        self._mass_engine_controls()
        self._mass_starter_pneumatic()
        self._mass_fuel_system()
        self._mass_flight_controls()
        self._mass_hydraulics()
        self._mass_electrical()
        self._mass_avionics()
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
        self.T  = convert.kgf_to_lbf(self.ac.propulsion.totalThrust)
        Vi = self.ac.designGoals.fuelMass/self.fuelProp['kerosene']
        self.Vi = convert.cubm_to_gal(Vi)
        
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
        We = convert.kg_to_lb(self.ac.propulsion.engine.mass)
        mEngineSection = 0.01*We**0.717*self.Ne*self.Nz
        return convert.lb_to_kg(mEngineSection), np.zeros(3)
        
    def _mass_airintake(self):
        #mAir = 13.29*Kvg * Ld**0.643* Kd**0.182 *Nen**1.498 * (Ls/Ld)**(-.373)*De
        return 0.0, np.zeros(3)
        
    def _mass_engine_cooling(self):
        #FIXME: engine shroud length is assumed 0.15*total length
        Lsh = convert.m_to_ft(self.ac.propulsion.engine.length)*0.15
        De = convert.m_to_ft(self.ac.propulsion.engine.diameter)
        mEc = 4.55*De*Lsh*self.Ne
        return convert.lb_to_kg(mEc), np.zeros(3)
        
    def _mass_oil_cooling(self):
        mOil = 37.82*self.Ne**1.023
        return convert.lb_to_kg(mOil), np.zeros(3)
    
    def _mass_engine_controls(self):
        #FIXME: assumed that engine controls are close to engine since it is UAV
        Lec = 0.25*convert.m_to_ft(self.ac.propulsion.engine.length)
        m = 10.5*self.Ne**1.008*Lec**0.222
        m = convert.lb_to_kg(m)
        self.total.airframe.add_item('Engine controls',m,np.zeros(3))
    
    def _mass_starter_pneumatic(self):
        Te = convert.kgf_to_lbf(self.ac.propulsion.engine.thrustMC)
        m = 0.025*Te**0.760*self.Ne**0.72
        m = convert.lb_to_kg(m)
        self.total.airframe.add_item('Engine starter',m,np.zeros(3))
    
    def _mass_fuel_system(self):
        #FIXME: assumption that all tanks are protected and pressurized
        Vi = 0.5*self.Vi
        Vt = self.Vi
        Vp = 0*self.Vi
        sfc = self.ac.propulsion.engine.sfcMC
        Nt = self.ac.propulsion.numberOfTanks
        m = 7.45*Vt**0.47*(1.+Vi/Vt)**(-0.095)*(1.+Vp/Vt)*Nt**0.066*self.Ne**0.052*(self.T*sfc/1e3)**0.249
        m = convert.lb_to_kg(m)
        self.total.airframe.add_item('Fuel system',m,np.zeros(3))
    
    def _mass_flight_controls(self):
        M = self.ac.designGoals.cruiseMach
        Scs = self.ac.wing.csArea
        Nc = 1.0 #FIXME: number of crew is set to 1 for UAV
        Ns = 2.0
        m = 36.28*M**0.003*Scs**0.489*Ns**0.484*Nc**0.127
        m = convert.lb_to_kg(m)
        self.total.airframe.add_item('Flight controls',m,np.zeros(3))
    
    def _mass_hydraulics(self):
        Kvsh = 1.0 #1.425 for variable sweep wing, otherwise 1.0
        Nu = 10.0 # number of hydraulic utility functions (typically 5-15)
        m = 37.23*Kvsh*Nu**0.664
        m = convert.lb_to_kg(m)
        self.total.airframe.add_item('Hydraulics',m,np.zeros(3))
    
    def _mass_electrical(self):
        Rkva = 120.0 # system electrical rating
        Kmc = 1.45
        Nc = 1.0
        La = convert.m_to_ft(self.ac.wing.chords[0])
        Ngen = self.Ne
        m = 172.2*Kmc*Rkva**0.152*Nc**0.1*La**0.1*Ngen**0.0091
        m = convert.lb_to_kg(m)
        self.total.airframe.add_item('Electronics',m,np.zeros(3))
    
    def _mass_avionics(self):
        mUav = 600.#FIXME: uninstalled avionics weight is assumed 1000lb
        m = 2.117*mUav**0.933
        m = convert.lb_to_kg(m)
        self.total.airframe.add_item('Avionics',m,np.zeros(3))

    