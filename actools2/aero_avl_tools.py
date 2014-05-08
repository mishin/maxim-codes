# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 12:49:55 2014

@author: Maxim
"""
def print_header(header,underline='-'):
    print str(header) + '\n'+underline*len(header)

from paths import MyPaths
import shlex
from subprocess import Popen, PIPE
import os

pth = MyPaths()


class AVL:
    def __init__(self):
        args = shlex.split(pth.avl,False,os.name=='posix')
        self.ps=Popen(args,stdin=PIPE,stderr=PIPE,stdout=PIPE)

    def cmd(self,command,echo=False):
        command = str(command)
        self.ps.stdin.write(command+'\n')
        if echo: print command
    def terminate(self):
        #self.cmd('\n\n\nQUIT')
        self.ps.stderr.close()
        self.ps.kill()

    def get_output(self):
        return str(self.ps.stdout.read()).replace('\r','')


class Results(object):
    def __init__(self):
        self.name = None
    def display(self):
        if self.name==None:
            print ''
        else:
            print_header(self.name)
        for attr, value in self.__dict__.iteritems():
            if type(value) is float:
                print '{:<8} = {:<+12.5f}  '.format(attr,value)


class AVLresults:
    def __init__(self,controlSurfaceIndex=['a','e','r','f']):
        self.csIndex = controlSurfaceIndex
        self.alpha=0.0
        self.beta=0.0
        self.Mach=0.0
        self.Sref=0.0
        self.Cref=0.0
        self.Bref=0.0
        self.ARref=0.0
        self.k=0.0
        self.e=0.0
        self.a=0.0
        self.xNP=0.0
        self.SM=0.0
        self.CD0=0.0
        self.coef = Coefficients()
        self.derivs = Derivatives(self.csIndex)
        self.hinges = HingeMoments(self.csIndex)
    def display(self):
        print_header('AVL analysis results','=')
        for attr, value in self.__dict__.iteritems():
            if type(value) is float:
                print '{:<8} = {:<+12.5f}  '.format(attr,value)
        self.coef.display()
        self.derivs.display()
        self.hinges.display()


class Derivatives(Results):
    def __init__(self,csIndex=['a','e','r','f']):
        """
        Parameters
        ----------
        
        csIndex : string array
            list of strings with control surface indexes. These index will be 
            used as control derivative name. For example a (aileron) -> CLda, CDda
            e (elevator) -> CLde, CDde
        """
        self.name = 'Derivatives'
        self.csIndex = csIndex

        self.CLa=0.0
        self.CYa=0.0
        self.Cla=0.0
        self.Cma=0.0
        self.Cna=0.0
        
        self.CLb=0.0
        self.CYb=0.0
        self.Clb=0.0
        self.Cmb=0.0
        self.Cnb=0.0
        
        self.CLp=0.0 
        self.CYp=0.0
        self.Clp=0.0
        self.Cmp=0.0
        self.Cnp=0.0
        
        self.CLp=0.0
        self.Clq=0.0
        self.CYq=0.0
        self.Cmq=0.0
        self.Cnq=0.0

        self.CLr=0.0
        self.CYr=0.0
        self.Clr=0.0
        self.Cmr=0.0
        self.Cnr=0.0
        self._init_cs_derivatives()
    
    def _init_cs_derivatives(self):
        for idx in self.csIndex:
            self.__dict__['CLd%s'%idx] = 0.0
            self.__dict__['CYd%s'%idx] = 0.0
            self.__dict__['Cld%s'%idx] = 0.0
            self.__dict__['Cmd%s'%idx] = 0.0
            self.__dict__['Cnd%s'%idx] = 0.0
            self.__dict__['CDd%s'%idx] = 0.0
            self.__dict__['ed%s'%idx]  = 0.0


class Coefficients(Results):
    def __init__(self):
        self.name = 'Coefficients'
        self.CX=0.0
        self.CY=0.0
        self.CZ=0.0
        self.Cl=0.0
        self.Cm=0.0
        self.Cn=0.0
        self.CL=0.0
        self.CD=0.0
        self.CDind=0.0


class HingeMoments(Results):
    def __init__(self,csIndex=['a','e','r','f']):
        self.name = 'Hinge moments'
        self.csIndex = csIndex
        self._init_cs()
    def _init_cs(self):
        for idx in self.csIndex:
            self.__dict__['H%s'%idx] = 0.0


def run_test1():
    rslt = AVLresults()
    rslt.display()

if __name__=="__main__":
    run_test1()