# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 13:28:36 2014

@author: Maxim
"""

import shlex
from subprocess import Popen, PIPE
import os
from airfoil_polar import AirfoilPolar
from paths import MyPaths
import numpy as np

pth = MyPaths()

def get_javafoil_analysis(airfoil,Mach,Re,alphaSeq):
    pth.set_file_prefix_random()
    tmpJournal = pth.get_tmp_file('jfscript')
    tmpAfFile = pth.get_tmp_file('dat')
    tmpPolar  = pth.get_tmp_file('pol')
    airfoil.write_txt(tmpAfFile)
    jouFile = open(tmpJournal,'wt')
    jouFile.write('Options.Country(0)\nGeometry.Clear()\n')
    jouFile.write('Geometry.Open(/%s/)\n' %tmpAfFile)
    jouFile.write('Options.MachNumber(%.4f)\n'%Mach)
    jouFile.write('Options.StallModel(0)\n')
    jouFile.write('Options.TransitionModel(1)\n')
    jouFile.write('Options.GroundEffect(0)\n')
    jouFile.write('Options.HeightOverSpan(0.5)\n')
    jouFile.write('Options.AspectRatio(0)\n')
    if os.name=='nt':
        jouFile.write('Polar.Analyze(%.0f;%.0f;%.0f;%.0f;%.0f;%.2f;100;100;0;0)\n'%(Re,Re,0,alphaSeq[0], alphaSeq[1], alphaSeq[2]))    
    else:
        jouFile.write('Polar.Analyze(%.0f:%.0f:%.0f:%.0f:%.0f:%.2f:100:100:0:0)\n'%(Re,Re,0,alphaSeq[0], alphaSeq[1], alphaSeq[2]))
    jouFile.write('Polar.Save(/%s/)\n'%tmpPolar)
    jouFile.write('Exit()')
    jouFile.close()
    cmd = ('%s -cp %s -jar %s Script=%s'%(pth.java,pth.mhclasses,pth.javafoil,tmpJournal))
    args=shlex.split(cmd,False,os.name=='posix')
    ps=Popen(args,stdin=PIPE,stderr=PIPE,stdout=PIPE)
    output, errors = ps.communicate()
    ps.stderr.close()
    polar = _read_javafoil_polar(tmpPolar)
    polar.Mach  = Mach
    polar.Re = Re
    os.remove(tmpJournal)
    os.remove(tmpAfFile)
    os.remove(tmpPolar)
    return polar


def _read_javafoil_polar(polarPath):
    """
    reads polar file of javafoil and returns AirfoilPolar object
    """
    polar = AirfoilPolar()
    polarFile = open(polarPath,'rt')
    lines = polarFile.readlines()
    polarFile.close()
    del lines[0:5]
    for line in lines:
        if line.strip()!='':
            segLine = line.split()
            try:
                polar.alpha = np.append(polar.alpha, float(segLine[0]))
                polar.cl    = np.append(polar.cl, float(segLine[1]))
                polar.cd    = np.append(polar.cd, float(segLine[2]))
                polar.cm    = np.append(polar.cm, float(segLine[3]))
                polar.TU    = np.append(polar.TU, float(segLine[4]))
                polar.TL    = np.append(polar.TL, float(segLine[5]))
                polar.SU    = np.append(polar.SU, float(segLine[6]))
                polar.SL    = np.append(polar.SL, float(segLine[7]))
                polar.LD    = np.append(polar.LD, float(segLine[8]))
                polar.AC    = np.append(polar.AC, float(segLine[9]))
                polar.CP    = np.append(polar.CP, float(segLine[10]))
            except ValueError:
                print 'Polar processing failed'
    return polar