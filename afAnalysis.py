# -*- coding: utf-8 -*-
"""
Created on Fri May 25 14:33:22 2012

@author: Maxim
"""
import subprocess as sp
import os
import afLib
import numpy as ny
from paths import myPaths
import scipy.optimize as root
import scipy.interpolate as interp
import matplotlib.pyplot as plt

class polar:
    def __init__(self,airfoil):
        self.airfoil = airfoil
        self.M  = []
        self.Re = []
        self.alpha = []
        self.CL = []
        self.CD = []
        self.CM = []
        self.TU = []
        self.TL = []

    def readXpolar(self,polarPath):
        class tmpPol:
            def __init__(self):
                self.alpha = []
                self.CL = []
                self.CD = []
                self.CDp = []
                self.CM = []
                self.TU = []
                self.TL = []

        tmpPolar = tmpPol()
        polarFile = open(polarPath,'rt')
        for ii in range(12): tmp = polarFile.readline()
        del tmp
        lines = polarFile.readlines()
        self.CDp = []
        for line in lines:
            segLine = line.split()
            tmpPolar.alpha = ny.append(tmpPolar.alpha,float(segLine[0]))
            tmpPolar.CL = ny.append(tmpPolar.CL,float(segLine[1]))
            tmpPolar.CD = ny.append(tmpPolar.CD,float(segLine[2]))
            tmpPolar.CDp = ny.append(tmpPolar.CDp,float(segLine[3]))
            tmpPolar.CM = ny.append(tmpPolar.CM,float(segLine[4]))
            tmpPolar.TU = ny.append(tmpPolar.TU,float(segLine[5]))
            tmpPolar.TL = ny.append(tmpPolar.TL,float(segLine[6]))
        polarFile.close()

        Sort = ny.argsort(tmpPolar.alpha)
        for ii in Sort:
            self.alpha = ny.append(self.alpha,tmpPolar.alpha[ii])
            self.CL = ny.append(self.CL,tmpPolar.CL[ii])
            self.CD = ny.append(self.CD,tmpPolar.CD[ii])
            self.CDp = ny.append(self.CDp,tmpPolar.CDp[ii])
            self.CM = ny.append(self.CM,tmpPolar.CM[ii])
            self.TU = ny.append(self.TU,tmpPolar.TU[ii])
            self.TL = ny.append(self.TL,tmpPolar.TL[ii])
        del tmpPolar
    
        if self.CL == []:
            self.CLmax = -100.
        else:
            self.CLmax = ny.max(self.CL)
            
    def readJpolar(self,filePath):
        polarFile = open(filePath,'rt')
        for ii in range(5): tmp = polarFile.readline()
        del tmp
        lines = polarFile.readlines()
        for ii in range(2): lines = ny.delete(lines,-1,0)
        polarFile.close()
        self.SU = []
        self.SL = []
        self.LD = []
        self.AC = []
        self.CP = []
        for line in lines:
            segLine = line.split()
            self.alpha = ny.append(self.alpha, float(segLine[0]))
            self.CL = ny.append(self.CL, float(segLine[1]))
            self.CD = ny.append(self.CD, float(segLine[2]))
            self.CM = ny.append(self.CM, float(segLine[3]))
            self.TU = ny.append(self.TU, float(segLine[4]))
            self.TL = ny.append(self.TL, float(segLine[5]))
            self.SU = ny.append(self.SU, float(segLine[6]))
            self.SL = ny.append(self.SL, float(segLine[7]))
            self.LD = ny.append(self.LD, float(segLine[8]))
            self.AC = ny.append(self.AC, float(segLine[9]))
            self.CP = ny.append(self.CP, float(segLine[10]))
        self.CLmax = ny.max(self.CL)
        self.ClmaxAlpha = self.alpha[ny.argmax(self.CL)]

    def calcXpolar(self,Mach, Re,alphaStart, alphaEnd, alphaStep, Iter=10, Graphic = True, Smooth=False):
        path = myPaths()
        path.setRandPrefix()
        tmpAfFile = path.getTmpFile('dat','af')
        tmpPolar  = path.getTmpFile('pol')
        tmpDump   = path.getTmpFile('dmp')

        afLib.writeAirfoil(tmpAfFile,self.airfoil)
        self.M = Mach
        self.Re = Re
    
        def issueCmd(cmd,echo=True):
            ps.stdin.write(cmd+'\n')
            if echo: print cmd
        ps = sp.Popen([path.Xfoil], stdin=sp.PIPE, stdout=None, stderr=None)
        
        if Graphic == False: issueCmd('PLOP\nG\n')
        issueCmd('LOAD')
        issueCmd('%s' %tmpAfFile)
        if Smooth == True:issueCmd('GDES\nCADD\n\n\n\n\nPANEL')
        issueCmd('OPER')
        issueCmd('VISC\n%.0f' % self.Re)
        issueCmd('MACH\n%.4f' % self.M)
        if Iter>10: issueCmd('ITER\n%d' % Iter)
        issueCmd('PACC')
        issueCmd('%s\n%s' % (tmpPolar,tmpDump))

        if alphaStep == 0:
            issueCmd('ALFA\n%.2f'%alphaStart)
        elif alphaStart*alphaEnd<0:
            issueCmd('ASEQ\n%.2f\n%.2f\n%.2f'%(-alphaStep, alphaStart, -alphaStep))
            issueCmd('ASEQ\n%.2f\n%.2f\n%.2f'%(0, alphaEnd, alphaStep))
        else:
            issueCmd('ASEQ\n%.2f\n%.2f\n%.2f'%(alphaStart, alphaEnd, alphaStep))
        issueCmd('\nQUIT')
        ps.wait()
    
        self.readXpolar(tmpPolar)
    
        os.remove(tmpAfFile)
        os.remove(tmpPolar)
        os.remove(tmpDump)
        
    def calcJpolar(self, Mach, Re, alphaStart, alphaEnd, alphaStep,flap=False):
        path = myPaths()
        path.setRandPrefix()
        tmpJournal = path.getTmpFile('jfscript')
        tmpAfFile = path.getTmpFile('dat','af')
        tmpPolar  = path.getTmpFile('pol')
        
        if flap:
            afLib.writeFlap(tmpAfFile, self.airfoil)
        else:
            afLib.writeAirfoil(tmpAfFile,self.airfoil)
            
        self.M = Mach
        self.Re = Re

        jouFile = open(tmpJournal,'wt')
        jouFile.write('Options.Country(0)\nGeometry.Clear()\n')
        jouFile.write('Geometry.Open(\"%s\")\n' %tmpAfFile)
        jouFile.write('Options.MachNumber(%.4f)\n'%Mach)
        jouFile.write('Options.StallModel(0)\n')
        jouFile.write('Options.TransitionModel(1)\n')
        jouFile.write('Options.GroundEffect(0)\n')
        jouFile.write('Options.HeightOverSpan(0.5)\n')
        jouFile.write('Options.AspectRatio(0)\n')
        jouFile.write('Polar.Analyze(%.0f;%.0f;%.0f;%.0f;%.0f;%.0f;100;100;0;0)\n'%(Re,Re,0,alphaStart, alphaEnd, alphaStep))
        jouFile.write('Polar.Save(\"%s\")\n'%tmpPolar)
        jouFile.write('Exit()')
        jouFile.close()
        cmd = ('\"\"%s\" -cp \"%s\" -jar \"%s\" Script=\"%s\"\"'%(path.java,path.mhclasses,path.javafoil,tmpJournal))
        os.system(cmd)
        self.readJpolar(tmpPolar)
        os.remove(tmpJournal)
        os.remove(tmpAfFile)
        os.remove(tmpPolar)

    def cdAtcl(self,clreq):
        try:
            clalpha = interp.interp1d(self.alpha,self.CL,'linear')
            cdalpha = interp.interp1d(self.alpha,self.CD,'cubic')
            def fcn(x): 
                x = float(x)
                return clalpha(x)-clreq
            alpha = root.bisect(fcn,-5.0,10.0)
            print 'alphareq', alpha
            cdreq = cdalpha(alpha)
        except:
            cdreq = float(100)
            alpha = 100
            print 'Drag coefficient', cdreq
        return cdreq, alpha

def testFcn():
    af = 'GA37A315mod.dat'
    airfoil = afLib.readAirfoil(af)
    Polar = polar(airfoil)
    Polar.calcJpolar(0.16,4e6,-30,30,1)
    
    plt.plot(Polar.alpha, Polar.CD)
    plt.show()

#testFcn()