# -*- coding: utf-8 -*-
"""
Created on Tue Nov 05 18:00:51 2013

@author: Maxim
"""
from numpy import array, zeros, linspace

class Mesh():
    def __init__(self):
        self._xfarfield = 20.0
        self._yfarfield = 15.0
        self._airfoilPts = 50
        self._xfarfieldPts = 50
        self._yfarfieldPts = 50
        self._wallSpacing = 1e-4
        self._farfieldSpacing = 1.0
        self._leSpacing = 1e-3
        self._teSpacing = 1e-3
    
    def create(self,meshFilePath):
        self._create_script()
        #run_scripts
    
    def _create_script(self,igsAirfoil=True):
        if igsAirfoil:
            open_igs()
        else:
            create_akima_airfoil()
        open_template()
        copy_lines_to_new_file()
        write_custom_line()
        close_file()