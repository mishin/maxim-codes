# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 11:59:44 2014

@author: Maxim
"""

from subprocess import Popen, PIPE
import shlex
import os
from airfoil_polar import AirfoilPolar, AirfoilPolar1D
from paths import MyPaths
import numpy as np

pth = MyPaths()

class Xfoil:
    def __init__(self,graphic=False):
        args = shlex.split(pth.Xfoil,False,os.name=='posix')
        self.ps=Popen(args,stdin=PIPE,stderr=PIPE,stdout=PIPE)
        if graphic==False:
            self.cmd('PLOP\nG\n')
    def cmd(self,command,echo=False):
        self.ps.stdin.write(command+'\n')
        if echo: print command
    def terminate(self):
        self.cmd('\n\n\nQUIT')
        self.ps.stderr.close()
        self.ps.stdout.read()


def get_xfoil_analysis(airfoil,Mach,Re,alphaSeq=[-10,20,1.0],
                       nIter=10,graphic=False,smooth=False):

    pth.set_file_prefix_random()
    tmpAfFile = pth.get_tmp_file('txt','_crd')
    tmpPolar  = pth.get_tmp_file('txt','_pol')
    airfoil.write_txt(tmpAfFile,False)
    xfoil = Xfoil()
    xfoil.cmd('LOAD\n%s'%tmpAfFile)
    if smooth:
        xfoil.cmd('GDES\nCADD\n\n\n\n\nPANEL')
    xfoil.cmd('OPER\nVISC\n%.0f\nMACH\n%.4f'%(Re,Mach))
    if nIter>10:
        xfoil.cmd('ITER\n%d' % nIter)
    xfoil.cmd('PACC')
    xfoil.cmd(' ')
    xfoil.cmd(' ')
    if alphaSeq[2] == 0 or alphaSeq[0]==alphaSeq[1]:
        xfoil.cmd('ALFA\n%.2f'%alphaSeq[0])
    elif alphaSeq[0]*alphaSeq[1]<0:
        xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(0, alphaSeq[1], alphaSeq[2]))
        xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(-alphaSeq[2], alphaSeq[0], -alphaSeq[2]))
    elif alphaSeq[0]*alphaSeq[1]>=0 and alphaSeq[0]>=0:
        xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(alphaSeq[0], alphaSeq[1], alphaSeq[2]))
    else:
        xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(alphaSeq[1], alphaSeq[0], -alphaSeq[2]))
    xfoil.cmd('PWRT')
    xfoil.cmd('%s'%(tmpPolar))
    xfoil.terminate()
    polar = _read_xfoil_polar(tmpPolar)
    polar.Mach = Mach
    polar.Re   = Re
    polar._calc_clmax()
    os.remove(tmpAfFile)
    os.remove(tmpPolar)
    return polar


def _read_xfoil_polar(polarPath):
    """
    reads polar file of Xfoil and returns AirfoilPolar object
    """
    polar = AirfoilPolar1D()
    polarFile = open(polarPath,'rt')
    lines = polarFile.readlines()
    polarFile.close()
    del lines[0:12]
    newLines = [float(s) for s in lines[0].split()]
    for line in lines[1::]:
        addLine = [float(s) for s in line.split()]
        newLines = np.vstack([newLines,addLine])
    newIdx      = np.argsort(newLines[:,0])
    polar.alpha = np.array(newLines[newIdx,0])
    polar.cl    = np.array(newLines[newIdx,1])
    polar.cd    = np.array(newLines[newIdx,2])
    polar.cdp  = np.array(newLines[newIdx,3])
    polar.cm    = np.array(newLines[newIdx,4])
    return polar


def create_by_interp_txt(afPath1,afPath2,afName='interpolated',
                         interpFraction=0.5,graphic=False):
    """
    Interpolates or "blending" airfoil in various proportions using Xfoil.
    
    Parameters
    ----------
    
    afPath1 : path
        file path of first airfoil text fle
    afPath2 : path
        file path of second airfoil text file
    afName : string
        name of result airfoil
    interpFraction : float [0,1]
        defines in which proportion airfoils will be merged. 0 - will 
        return 1st airfoil, 1 - will return 2nd airfoil.
    graphic : bool
        enable/disable Xfoil gui appearance
    """
    xfoil = Xfoil()
    pth.set_file_prefix_random()
    tmpAfFile = pth.get_tmp_file('dat')
    xfoil.cmd('INTE')
    xfoil.cmd('F')
    xfoil.cmd('%s'%afPath1)
    xfoil.cmd('F')
    xfoil.cmd('%s'%afPath2)
    xfoil.cmd('%.4f'%interpFraction)
    xfoil.cmd('%s'%afName)
    xfoil.cmd('PCOP')
    xfoil.cmd('SAVE')
    xfoil.cmd('%s'%tmpAfFile)
    xfoil.terminate()
    return tmpAfFile


if __name__=="__main__":
    print 'hello world!'