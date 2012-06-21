# -*- coding: utf-8 -*-
"""
Created on Mon Jun 04 10:25:42 2012

@author: Maxim
"""
import os

class myPaths:
    def __init__(self):
        self.wdir      = os.getcwd()
        self.Xfoil     = self.wdir + '\Xfoil\Xfoil.exe'
        self.javafoil  = self.wdir + '\Jfoil\javafoil.jar'
        self.mhclasses = self.wdir +'\Jfoil\mhclasses.jar'
        self.tmpdir    = self.wdir + '\\temp'
        self.db        = self.wdir + '\\db'
        self.java      = r'C:\Program Files\Java\jre7\bin\java.exe'

    def setRandPrefix(self):
        from random import randrange
        self.NamePrefix = str(randrange(0,1e6,1))

    def getTmpFile(self,fileExt = "", addSymbol = ""):
        fileName = self.tmpdir + '\\' + self.NamePrefix + addSymbol
        if fileExt != "" : fileName  = fileName + "." + fileExt
        return fileName