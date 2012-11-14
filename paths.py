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

    def setRandPrefix(self,N=7):
        from random import randrange
        prefix = ''
        for ii in range(N):
            r = randrange(0,3)
            if r == 0:
                prefix += str(randrange(0,10))
            elif r == 1:
                prefix += chr(randrange(65,91))
            else:
                prefix += chr(randrange(97,123))
        self.NamePrefix = prefix

    def getTmpFile(self,fileExt = "", addSymbol = ""):
        fileName = self.tmpdir + '\\' + self.NamePrefix + addSymbol
        if fileExt != "" : fileName  = fileName + "." + fileExt
        return fileName
        
def testFcn():
    a = myPaths()
    for ii in range(5):
        a.setRandPrefix()
        print ii, '\t',a.NamePrefix,'\n'
    
#testFcn()