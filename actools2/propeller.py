# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 17:11:11 2014

@author: Maxim
"""
from paths import MyPaths
import db_tools
import airfoil

pth = MyPaths()

def load(propellerName):
    pass

class Propeller:
    def __init__(self):
        self.name = ''
        self.diameter     = 0.0
        self.radius       = 0.0
        self.diameterHub  = 0.0
        self.numBlades    = 3
        self._numSections = 0
        self.r            = list()
        self.x            = list()
        self.chord        = list()
        self.beta         = list()
        self.airfoilNames = list()
        self.airfoils     = list()

    def read_db(self,name,xlsPath=None,airfoilXlsPath=None):
        if xlsPath==None:
            xlsPath = pth.db.prop
        if airfoilXlsPath==None:
            airfoilXlsPath = pth.db.propAirfoil
        sh = db_tools.ReadDatabase(xlsPath,name)
        self.name = name
        self.diameter     = sh.read_row(1,1)
        self.diameterHub  = sh.read_row(-1,1)
        self.numBlades    = sh.read_row(-1,1)
        self.betaRange    = sh.read_row(-1,1)
        self.r            = sh.read_row(-1,1)
        self.chord        = sh.read_row(-1,1)
        self.beta         = sh.read_row(-1,1)
        self.airfoilNames = sh.read_row(-1,1)
        self._numSections = len(self.r)
        for afname in self.airfoilNames:
            self.airfoils.append(airfoil.load(afname, airfoilXlsPath))
        self._analyze_geometry()

    def write_db(self,xlsPath=None):
        if xlsPath==None:
            xlsPath = pth.db.prop
    
    def _create_splines(self):
        pass
    
    def _analyze_geometry(self):
        self.radius = self.diameter / 2.0
        self.x = self.r / self.radius
        self._calc_solidity_ratio()
        
    def set_beta(self,betaNew):
        pass
    
    def analyze(self,beta,rpm,velocity,density):
        pass
    
    def display(self):
        pass