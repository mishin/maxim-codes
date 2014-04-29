# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 21:39:57 2014

@author: Maxim
"""
import numpy as np
import constants
from weight_tools import AircraftMass, AircraftMass2
import convert


def get_flying_wing_mass(aircraft):
    mass = BlendedWingBodyMass(aircraft)
    mass.analyze()
    mass.total.fuel = aircraft.mass.fuel
    mass.total.payload = aircraft.mass.payload
    mass.total.update_total()
    mass.output.display()
    return mass.total


class BlendedWingBodyMass(object):
    def __init__(self,aircraft):
        self.ac = aircraft
        self.output = AircraftMass2(self.ac.name,self.ac.wing.MAC,self.ac.wing.MAClocation[0])
        self.total = AircraftMass(self.ac.name,self.ac.wing.MAC,self.ac.wing.MAClocation[0])
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
        #self.total.display()
    
    def _get_coefficients(self):
        self.Kdw  = 1.0 # 0.786 for delta wing, otherwise Kdw=1.0
        self.Kvs  = 1.0 # 1.19 for variable sweep, otherwise 1.0
        self.Kdwf = 1.0 # 0.774 for delta wing, otherwise 1.0
        self.Kcb  = 1.0
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
    
    def _add_mass1(self,name,mass,cg=None):
        """ NOTE: CG is in meters"""
        m = convert.lb_to_kg(mass)
        self.total.airframe.add_item(name, m, cg)
        self.output.empty.add_item(name,m,cg)

    def _mass_wing(self,wing):
        # by Raymer - fighter weights
        Sw = convert.sqm_to_sqft(wing.area)
        A  = wing.aspectRatio
        tcRoot = wing.airfoils[0].thickness
        lmbda = wing.sweepElasticRad #FIXME: check rad/deg
        Scsw  = 0.1*Sw #FIXME: calculate control surface area
        m = 0.0103*self.Kdw*self.Kvs* (self.Wdg*self.Nz)**0.5* Sw**0.622* A**0.785 *tcRoot**(-0.4)* (1+lmbda)**0.5*np.cos(lmbda)**(-1.0)* Scsw**0.04
        wingCGratio = self.constMass['wingCGratio']
        xCG = self.ac.wing.MAClocation[0] + self.ac.wing.MAC*wingCGratio
        zCG = self.ac.wing.secApex[0,2]
        self._add_mass1('wing',m,np.array([xCG,0.0,zCG]))
    
    def _mass_fuselage(self):
        L = convert.m_to_ft(self.ac.wing.chords[0])
        D = convert.m_to_ft(self.ac.wing.airfoils[0].thickness)
        W = convert.m_to_ft(self.ac.fusWidth)
        m = 0.499*self.Kdwf*self.Wdg**0.35 *self.Nz**0.25 *L**0.5 *D**0.849 *W**0.685
        xCG = self.ac.wing.chords[0] * self.constMass['fuseCGratio'][1]
        zCG = self.ac.wing.secApex[0,2]
        self._add_mass1('fuselage',m,np.array([xCG,0.0,zCG]))
    
    def _mass_mlg(self):
        Lm = convert.m_to_ft(self.ac.landingGear.strutLength[1]) #TODO: check value
        m = self.Kcb*self.Ktpg*(self.Wt*self.Nl)**0.25 * Lm**0.973
        xCG = self.ac.landingGear.groundContactX[1]
        self._add_mass1('main landing gear',m,np.array([xCG,0,0]))

    def _mass_nlg(self):
        Ln = convert.m_to_ft(self.ac.landingGear.strutLength[0])
        Nnw = 1.0 # number of nose wheels
        m = (self.Wl*self.Nl)**0.29 *Ln**0.5 *Nnw**0.525
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
        #mAir = 13.29*Kvg * Ld**0.643* Kd**0.182 *Nen**1.498 * (Ls/Ld)**(-.373)*De
        return 0.0, np.zeros(3)
        
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
        mUav = 600.#FIXME: uninstalled avionics weight is assumed 1000lb
        m = 2.117*mUav**0.933
        self._add_mass1('avionics',m)

    