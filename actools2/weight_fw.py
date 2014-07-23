# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 21:39:57 2014

@author: Maxim
"""
import numpy as np
import constants
from weight_tools import AircraftMass, MassList, Fuel
import convert

def get_flying_wing_mass(aircraft):
    out = AircraftMass(aircraft.name,aircraft.wing.MAC,aircraft.wing.MAClocation[0])
    mass = BlendedWingBodyMass(aircraft)
    mass.analyze()
    out.empty = mass.emptyMass
    out.payload = aircraft.mass.payload
    out.fuel = aircraft.mass.fuel
    out.update_total()
    return out


class BlendedWingBodyMass(object):
    def __init__(self,aircraft):
        self.ac = aircraft
        self.emptyMass = MassList(self.ac.name)
        self.fuelProp  = constants.load('fuel_density')
        self.constMass = constants.load('mass')
        
    def analyze(self):
        self._get_data()
        self._mass_wing(self.ac.wing)
        self._mass_fuselage()
        self._mass_mlg()
        self._mass_nlg()
        self._mass_engine_mount()
        self._mass_engine_section()
        self._mass_engine_cooling()
        self._mass_airintake()
        self._mass_oil_cooling()
        self._mass_engine_controls()
        self._mass_starter_pneumatic()
        self._mass_fuel_system()
        self._mass_flight_controls()
        self._mass_hydraulics()
        self._mass_electrical()
        self._mass_avionics()
        self._mass_engine()
    
    def _get_coefficients(self):
        self.Kdw  = 1.0 # 0.786 for delta wing, otherwise Kdw=1.0
        self.Kvs  = 1.0 # 1.19 for variable sweep, otherwise 1.0
        self.Kdwf = 1.0 # 0.774 for delta wing, otherwise 1.0
        self.Kcb  = 1.0 # 2.25 for cross beam landing gear, otherwise 1.0
        self.Ktpg = 1.0 # 0.826 for tripod landing gear, otherwise 1.0

    def _get_data(self):
        self._get_coefficients()
        self.Nz = 1.5 * self.ac.designGoals.loadFactor # ultimate load factor
        self.Wdg = convert.kg_to_lb(self.ac.designGoals.grossMass)
        self.Nl = self.ac.designGoals.loadFactorLanding
        self.Wt = self.Wdg
        self.Wl = self.Wdg
        self.Ne = self.ac.propulsion.numberOfEngines
        self.T  = convert.kgf_to_lbf(self.ac.propulsion.totalThrust)
        self.mUav = convert.kg_to_lb(self.ac.designGoals.avionicsMass)
        Vi = self.ac.designGoals.fuelMass/self.fuelProp['kerosene']
        self.Vi = convert.cubm_to_gal(Vi)
    
    def _add_mass1(self,name,mass,cg=None,lb=True):
        """ NOTE: CG is in meters"""
        if lb:
            mass = convert.lb_to_kg(mass)
        self.emptyMass.add_item(name, mass, cg)
    
    def _mass_engine(self):
        mKg = self.ac.propulsion.totalMass
        CG = self.ac.propulsion.CG
        self._add_mass1('engine',mKg,CG,False)

    def _mass_wing(self,wing):
        # by Raymer - fighter weights
        Sw = convert.sqm_to_sqft(wing.area)
        A  = wing.aspectRatio
        tcRoot = wing.airfoils[1].thickness #FIXME: assume thickness of 2nd section
        # --- own corrections ---
        eqArea = (wing.chords[0]+wing.chords[-1])*wing.span/2.0
        corr = wing.area/eqArea
        Kms = 1.1 # correction for multisegment wings
        TR = wing.taperRatio
        sweepC4 = wing.equivSweepC4rad
        Scsw  = wing.csArea
        try:
            m = 0.0103*self.Kdw*self.Kvs* (self.Wdg*self.Nz)**0.5* Sw**0.622* A**0.785 *tcRoot**(-0.4)* (1+TR)**0.5*np.cos(sweepC4)**(-1.0)* Scsw**0.04
        except ValueError:
            print Sw, A, tcRoot, TR, sweepC4, Scsw
            self.ac.display()
            raw_input()
        m = corr*Kms*m
        wingCGratio = self.constMass['wingCGratio']
        xCG = self.ac.wing.MAClocation[0] + self.ac.wing.MAC*wingCGratio
        zCG = self.ac.wing.secApex[0,2]
        if wing.material=='composite':
            m = m*self.constMass['compositeWing']
        self._add_mass1('wing',m,np.array([xCG,0.0,zCG]))
    
    def _mass_fuselage(self):
        L = convert.m_to_ft(self.ac.wing.chords[0])
        D = convert.m_to_ft(self.ac.wing.airfoils[0].thickness)
        W = convert.m_to_ft(self.ac.fusWidth)
        m = 0.499*self.Kdwf*self.Wdg**0.35 *self.Nz**0.25 *L**0.5 *D**0.849 *W**0.685
        xCG = self.ac.wing.chords[0] * self.constMass['fuseCGratio'][1]
        zCG = self.ac.wing.secApex[0,2]
        if self.ac.wing.material=='composite':
            m = m*self.constMass['compositeBody']
        self._add_mass1('fuselage',m,np.array([xCG,0.0,zCG]))
    
    def _mass_mlg(self):
#        Lm = convert.m_to_ft(self.ac.landingGear.strutLength[1]) #TODO: check value
#        m = self.Kcb*self.Ktpg*(self.Wt*self.Nl)**0.25 * Lm**0.973
        m = self.Wt * 0.04*0.85
        xCG = self.ac.landingGear.groundContactX[1]
        self._add_mass1('main landing gear',m,np.array([xCG,0,0]))

    def _mass_nlg(self):
#        Ln = convert.m_to_ft(self.ac.landingGear.strutLength[0])
#        Nnw = 1.0 # number of nose wheels
#        m = (self.Wl*self.Nl)**0.29 *Ln**0.5 *Nnw**0.525
        m = self.Wt * 0.04*0.15
        xCG = self.ac.landingGear.groundContactX[0]
        self._add_mass1('nose landing gear',m,np.array([xCG,0,0]))

    def _mass_engine_mount(self):
        m = 0.013*self.Ne**0.795*self.T**0.579 *self.Nz
        self._add_mass1('engine mount',m)
        
    def _mass_engine_section(self):
        We = convert.kg_to_lb(self.ac.propulsion.engine.mass)
        m = 0.01*We**0.717*self.Ne*self.Nz
        self._add_mass1('engine section',m)
        
    def _mass_airintake(self):
        Kvg = 1.0 # 1.62 for variable geometry
        Kd = 3.43 # duct constant. see fig 15.2 Raymer
        Ld = self.ac.wing.chords[0]-self.ac.propulsion.engine.length
        Ld = convert.m_to_ft(Ld)
        Ls = Ld
        De = convert.m_to_ft(self.ac.propulsion.engine.diameter)
        mAir = 13.29*Kvg * Ld**0.643* Kd**0.182 *self.Ne**1.498 * (Ls/Ld)**(-.373)*De
        self._add_mass1('air intake',mAir)
        
    def _mass_engine_cooling(self):
        #FIXME: engine shroud length is assumed 0.15*total length
        Lsh = convert.m_to_ft(self.ac.propulsion.engine.length)*0.15
        De = convert.m_to_ft(self.ac.propulsion.engine.diameter)
        m = 4.55*De*Lsh*self.Ne
        self._add_mass1('engine cooling',m)
        
    def _mass_oil_cooling(self):
        m = 37.82*self.Ne**1.023
        self._add_mass1('oil cooling',m)
    
    def _mass_engine_controls(self):
        #FIXME: assumed that engine controls are close to engine since it is UAV
        Lec = 0.25*convert.m_to_ft(self.ac.propulsion.engine.length)
        m = 10.5*self.Ne**1.008*Lec**0.222
        self._add_mass1('engine controls',m)
    
    def _mass_starter_pneumatic(self):
        Te = convert.kgf_to_lbf(self.ac.propulsion.engine.thrustMC)
        m = 0.025*Te**0.760*self.Ne**0.72
        self._add_mass1('engine starter',m)
    
    def _mass_fuel_system(self):
        #FIXME: assumption that all tanks are protected and pressurized
        Vi = 0.5*self.Vi
        Vt = self.Vi
        Vp = 0*self.Vi
        sfc = self.ac.propulsion.engine.sfcMC
        Nt = self.ac.propulsion.numberOfTanks
        m = 7.45*Vt**0.47*(1.+Vi/Vt)**(-0.095)*(1.+Vp/Vt)*Nt**0.066*self.Ne**0.052*(self.T*sfc/1e3)**0.249
        self._add_mass1('fuel system',m)
    
    def _mass_flight_controls(self):
        M = self.ac.designGoals.cruiseMach
        Scs = self.ac.wing.csArea
        Nc = 1.0 #FIXME: number of crew is set to 1 for UAV
        Ns = 2.0
        m = 36.28*M**0.003*Scs**0.489*Ns**0.484*Nc**0.127
        self._add_mass1('flight controls',m)
    
    def _mass_hydraulics(self):
        Kvsh = 1.0 #1.425 for variable sweep wing, otherwise 1.0
        Nu = 10.0 # number of hydraulic utility functions (typically 5-15)
        m = 37.23*Kvsh*Nu**0.664
        self._add_mass1('hydraulics',m)
    
    def _mass_electrical(self):
        Rkva = 120.0 # system electrical rating
        Kmc = 1.45
        Nc = 1.0
        La = convert.m_to_ft(self.ac.wing.chords[0])
        Ngen = self.Ne
        m = 172.2*Kmc*Rkva**0.152*Nc**0.1*La**0.1*Ngen**0.0091
        self._add_mass1('electric system',m)
    
    def _mass_avionics(self):
        #mUav = 600.#FIXME: uninstalled avionics weight is assumed 1000lb
        mUav = self.mUav
        m = 2.117*mUav**0.933
        self._add_mass1('avionics',m)

    