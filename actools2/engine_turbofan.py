# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:07:43 2014

@author: Maxim
"""
from db_tools import ReadDatabase
from flight_conditions import ISAtmosphere
from paths import MyPaths
from weight_tools import get_total_cg
from numpy import ones,array,linspace
from math import isnan
from engine_turbofan_analysis import engine_modeling
from scipy.interpolate import interp1d, Rbf
import matplotlib.pyplot as plt
import constants
from misc_tools import Normalization, read_tabulated_data_with_header
import os.path


class Propulsion(object):
    def __init__(self):
        self.numberOfEngines = 1
        self.engine = TurbofanEngine()
        self.CGx = None
        self.CGy = None
        self.CGz = None
        self.numberOfTanks = 1
        self.sfcModel = None

    def load(self,name,buildTable=False):
        self.engine.load(name)
        self._process_data(buildTable)

    def _process_data(self,buildTable=False):
        self.totalThrust = self.engine.thrustMC*self.numberOfEngines
        self._calc_cg()
        if buildTable:
            path = MyPaths()
            thrustTablePath = '%s//ThrustTable_%s.txt'%(path.db.wdir,self.engine.name)
            if os.path.isfile(thrustTablePath):
                #print 'loading thrust table'
                self._load_thrust_table(thrustTablePath)
            else:
                #print 'building thrust table'
                self._build_thrust_table()

    def _calc_cg(self):
        if hasattr(self.CGx,'__iter__'):
            mass = ones(self.numberOfEngines)*self.engine.mass
            self.totalMass, self.CG = get_total_cg(mass,self.CGx,self.CGy,self.CGz)


    def get_sfc_direct(self,Mach,altitude,thrustReq):
        Mach      = float(Mach)
        altitude  = float(altitude)
        thrustReq = float(thrustReq) / self.numberOfEngines # thrust per engine
        
        sfcPerEngine = self.engine.get_sfc(Mach,altitude,thrustReq)
        return sfcPerEngine*self.numberOfEngines /3600. # kg/(N*sec)
    
    def get_sfc2(self,Mach,altitude,thrustReq):
        Mach = self._normMach.normalize(Mach)
        alt = self._normAlt.normalize(altitude)
        Treq = self._normTreq.normalize(thrustReq)
        return self.sfcModel(Mach,alt,Treq)
    
    def get_sfc(self,Mach,altitude,thrustReq):
        if self.sfcModel==None:
            return self.get_sfc_direct(Mach,altitude,thrustReq)
        else:
            return self.get_sfc2(Mach,altitude,thrustReq)
    
    def _load_thrust_table(self,path):
        self._get_thrust_table_limits()
        data = read_tabulated_data_with_header(path)
        self._normMach = Normalization(data['Mach'].min(), data['Mach'].max(),0,1)
        self._normAlt = Normalization(data['altitude'].min(), data['altitude'].max(),0,1)
        self._normTreq = Normalization(data['ThrustRequired'].min(), data['ThrustRequired'].max(),0,1)
        normMach = self._normMach.normalize(data['Mach'])
        normAlt = self._normAlt.normalize(data['altitude'])
        normTreq = self._normTreq.normalize(data['ThrustRequired'])
        self.sfcModel = Rbf(normMach,normAlt,normTreq,data['sfc'])
    
    def _get_thrust_table_limits(self):
        n = 10
        MachMin = 0.05
        MachMax = 1.0
        altMin = 0
        altMax = 2e4
        Tmc = self.totalThrust
        Tmin = 0.05 *Tmc
        Tmax = Tmc
        self._normMach = Normalization(MachMin,MachMax, 0,1)
        self._normAlt  = Normalization(altMin,altMax,0,1)
        self._normTreq = Normalization(Tmin, Tmax,0,1)
        return n, MachMin, MachMax, altMin, altMax, Tmin, Tmax

    def _build_thrust_table(self):
        n, MachMin, MachMax, altMin, altMax, Tmin, Tmax = self._get_thrust_table_limits()
        MachList = linspace(MachMin, MachMax, n)
        altList  = linspace(altMin, altMax, n)
        TreqList = linspace(Tmin, Tmax, n)
        nsfc = n*n*n
        sfc       = ones(nsfc)
        Mach      = ones(nsfc)
        altitude  = ones(nsfc)
        thrustReq = ones(nsfc)
        
        i = 0
        for M in MachList:
            for alt in altList:
                for Treq in TreqList:
                    _sfc = self.get_sfc_direct(M,alt,Treq)
                    if not isnan(_sfc):
                        sfc[i] = _sfc
                        Mach[i] = M #self._normMach.normalize(M)
                        altitude[i] = alt #self._normAlt.normalize(alt)
                        thrustReq[i] = Treq #self._normTreq.normalize(Treq)
                        i += 1
                        #print '%d\t%.4f\t%.0f\t%.0f\t%.6f'%(i, M, alt, Treq, _sfc)
        
        path = MyPaths()
        thrustTablePath = '%s//ThrustTable_%s.txt'%(path.db.wdir,self.engine.name)
        fid = open(thrustTablePath,'wt')
        fid.write('Mach\taltitude\tThrustRequired\tsfc\n')
        for M,alt,Treq,_sfc in zip(Mach[:i],altitude[:i],thrustReq[:i],sfc[:i]):
            fid.write('%.32f\t%.32e\t%.32e\t%.32e\n'%(M,alt,Treq,_sfc))
        fid.close()

        normMach = self._normMach.normalize(Mach[:i])
        normAlt  = self._normAlt.normalize(altitude[:i])
        normTreq = self._normTreq.normalize(thrustReq[:i])
        self.sfcModel = Rbf(normMach, normAlt, normTreq, sfc[:i])

    def __call__(self,Mach,altitude,powerSetting):
        """ NOTE: powerSetting is linear function - will be replaced """
        thrustReq = self.engine.thrustMC * float(powerSetting)/100.
        return self.get_sfc(Mach,altitude,thrustReq)
    
    def get_thrust_available(self,altitude):
        atm = ISAtmosphere(altitude)
        rho0 = 1.2255
        return atm.density/rho0 *self.totalThrust

    def test(self):
        plt.figure()
        plt.hold(True)
        legend = list()
        alt = 1e4
        Mach = linspace(0.1,1,200)
        Tmc = self.totalThrust
        Pset = linspace(0.05*Tmc, Tmc ,5)
        sfc = ones([len(Pset),len(Mach)])
        sfc2 = ones([len(Pset),len(Mach)])

        for i,p in enumerate(Pset):
            for j,M in enumerate(Mach):
                sfc[i,j] = self.get_sfc_direct(M,alt,p)
                sfc2[i,j] = self.get_sfc2(M,alt,p)

            plt.plot(Mach,sfc[i],marker='*')
            plt.plot(Mach,sfc2[i])
            
            legend.append('%.0f'%p)
            legend.append('%.0frbf'%p)
        plt.legend(legend)
        plt.show()


class TurbofanEngine(object):
    def __init__(self):
        self.name         = None
        self.thrustMC     = 0.0
        self.sfcMC        = 0.0
        self.mass         = 0.0
        self.length       = 0.0
        self.diameter     = 0.0
        # for thrust analysis
        self.designMach    = 0.0
        self.deignAltitude = 0.0
        self.designThrust  = 0.0
        T = array([362.873896304, 2267.9618519, 22679.618519])
        ratio = array([3.4,5.5,8])
        self._massCurve = interp1d(T,ratio,'linear')
        # for thrust table
        
    
    def load(self,name, xlsPath=None):
        """ loads turbofan engine from database file"""
        path = MyPaths()
        if xlsPath==None:
            xlsPath = path.db.engineTurbofan
        self.name = str(name)
        db = ReadDatabase(xlsPath,self.name)
        self.thrustMC       = db.read_row(-1,1) *constants.GRAVITY_ACCEL
        self.sfcMC          = db.read_row(-1,1)
        self.length         = db.read_row(-1,1)
        self.diameter       = db.read_row(-1,1)
        self.mass           = db.read_row(-1,1,False)
        self._calc_data()
    
    def _calc_data(self):
        if self.mass==0.0:
            self.mass = self._massCurve(self.thrustMC)
        if self.designThrust==0.0:
            self.designThrust = self.thrustMC
    
    def display(self):
        print self.name
        print self.thrustMC
        print self.mass
    
    def get_sfc(self,Mach,altitude,thrustReq):
        Tstatic = self.thrustMC
        Mcr = self.designMach
        altD = self.deignAltitude
        d, Td, sfcD, sfcF = engine_modeling(Tstatic,Mcr,altD,thrustReq,altitude,Mach)
        return sfcF


def run_test3():
    prop = Propulsion()
    prop.load('F404',True)
    prop.test()

if __name__=="__main__":
    run_test3()
        