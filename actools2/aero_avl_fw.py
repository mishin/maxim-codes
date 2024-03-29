# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 13:06:22 2014

Main class for flying wing configuration aerodynamic analysis using AVL code.
This class can be used as parent for different configurations analysis.

1. code for flying wing analysis with elevator only

@author: Maxim
"""

from aero_avl_tools import AVLresults, AVL
from paths import MyPaths
import numpy as np
from flight_conditions import FlightConditions
import re

pth = MyPaths()

class FlightConditionsAVL(FlightConditions):
    def __init__(self,aircraft,velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                 mass=None,cg=None,inertia=None,CD0=None):
        self.Lref = aircraft.wing.MAC
        self.bref = aircraft.wing.span
        self.Sref = aircraft.wing.area
        if mass==None:
            mass = aircraft.get_mass()
        if cg==None:
            cg = aircraft.get_cg()
        if inertia==None:
            inertia =aircraft.get_inertia()
        if CD0==None:
            CD0 = aircraft.get_drag(velocity,altitude)
        self.set_flight_conditions(velocity,altitude,CmTrim,loadFactor,mass,cg,inertia,CD0)
    
    def set_flight_conditions(self,velocity,altitude,CmTrim=0.0,loadFactor=1.0,
                 mass=None,cg=None,inertia=None,CD0=None):
        Lref = self.Lref
        velocity = float(velocity)
        altitude = float(altitude)
        super(FlightConditionsAVL,self).__init__(velocity,altitude,0.0,Lref)
        self.CmTrim     = CmTrim
        self.loadFactor = loadFactor
        if not mass==None:
            self.mass = mass
        if not cg==None:
            self.cg = cg
        if not inertia==None:
            self.inertia = inertia
        if not CD0==None:
            self.CD0 = CD0

    def get_CLreq(self):
        W = self.mass * self.loadFactor * self.g
        rho = self.atm.density
        V = self.velocity
        return W/ (rho*V*V*self.Sref)


class AVLsolver(object):
    def __init__(self,aircraft):
        self.ac = aircraft

    def create_input_file(self):
        """
        basic terms
        wing configuration
        control surface - ELEVATOR only for now
        weight and CG
        """
        #pth.set_file_prefix('flyingwing')
        pth.set_file_prefix_random()
        path = pth.get_tmp_file('avl')
        cg = self.ac.get_cg()
        CD0 = self.ac.get_drag()
        fid = open(path,'wt')
        # general input
        fid.write('%s\n'%self.ac.name)
        fid.write('0.0\n')
        fid.write('0  0  0.0\n')
        fid.write('%.6f  %.6f  %.6f\n'%(self.ac.wing.area,self.ac.wing.MAC,self.ac.wing.span))
        fid.write('%.6f  %.6f  %.6f\n'%(cg[0],cg[1],cg[2]))
        fid.write('%.8f\n'%CD0)
        # wing input
        fid.write('SURFACE\n')
        fid.write('WING\n')
        fid.write('%d  0  %d  0\n'%(self.ac.vlm.panelsChordwise, self.ac.vlm.panelsSpanwise))
        fid.write('YDUPLICATE\n0\n')
        for i in range(self.ac.wing.nSec):
            fid.write('#--- section %d ---\n'%(i+1))
            pathAf = pth.get_tmp_file('txt','af%d'%(i+1))
            self.ac.wing.airfoils[i].write_txt(pathAf)
            apex = self.ac.wing.secApex[i]
            chord = self.ac.wing.chords[i]
            angle = self.ac.wing.secAngles[i]
            fid.write('SECTION\n')
            fid.write('%.12f  %.12f  %.12f  %.12f  %.12f\n'%(apex[0],apex[1],apex[2],chord,angle))
            fid.write('AFIL\n%s\n'%pathAf)
            fid.write('CONTROL\n')
            fid.write('elevator  1.0  %.8f  0  0  0  +1\n'%(self.ac.wing.elevon.location[i]))
        fid.close()
    
    def run_trim(self,fc):
        self.create_input_file()
        self.avl = AVL()
        self.run_avl(fc)
        result = self._process_output()
        result.CD0 = float(fc.CD0)
        result.SM = float((result.xNP - fc.cg[0])/fc.Lref)
        self.avl.terminate()
        #pth.clean_temp_files()
        return result

    def run_single_point(self,fc,alpha=0,beta=0,elevator=0):
        alpha    = float(alpha)
        beta     = float(beta)
        elevator = float(elevator)
        self.create_input_file()
        self.avl = AVL()
        self.run_avl(fc,False,alpha,beta,elevator)
        result = self._process_output()
        result.CD0 = float(fc.CD0)
        result.SM = float((result.xNP - fc.cg[0])/fc.Lref)
        self.avl.terminate()
        #pth.clean_temp_files()
        return result

    def run_avl(self,fc,runTrim=True,alpha=0,beta=0,elevator=0):
        self.avl.cmd('LOAD\n%s'%pth.get_tmp_file('avl'))
        self.avl.cmd('OPER')
        self.avl.cmd('O')
        self.avl.cmd('P')
        self.avl.cmd('T,T,T,T')
        self.avl.cmd('H')
        self.avl.cmd('D')
        self.avl.cmd(' ')
        self.avl.cmd('M')
        self.avl.cmd('M')
        self.avl.cmd(fc.mass*fc.loadFactor)
        self.avl.cmd('MN')
        self.avl.cmd(fc.Mach)
        self.avl.cmd('V')
        self.avl.cmd(fc.velocity)
        self.avl.cmd('D')
        self.avl.cmd(fc.atm.density)
        self.avl.cmd('G')
        self.avl.cmd(fc.g)
        self.avl.cmd('CD')
        self.avl.cmd(fc.CD0)
        self.avl.cmd('IX')
        self.avl.cmd(fc.inertia[0])
        self.avl.cmd('IY')
        self.avl.cmd(fc.inertia[1])
        self.avl.cmd('IZ')
        self.avl.cmd(fc.inertia[2])
        self.avl.cmd('X')
        self.avl.cmd(fc.cg[0])
        self.avl.cmd('Y')
        self.avl.cmd(fc.cg[1])
        self.avl.cmd('Z')
        self.avl.cmd(fc.cg[2])
        self.avl.cmd(' ')
        self.avl.cmd('C1')
        self.avl.cmd('V')
        self.avl.cmd(fc.velocity)
        self.avl.cmd('M')
        self.avl.cmd(fc.mass*fc.loadFactor)
        self.avl.cmd('D')
        self.avl.cmd(fc.atm.density)
        self.avl.cmd('G')
        self.avl.cmd(fc.g)
        self.avl.cmd('X')
        self.avl.cmd(fc.cg[0])
        self.avl.cmd('Y')
        self.avl.cmd(fc.cg[1])
        self.avl.cmd('Z')
        self.avl.cmd(fc.cg[2]) 
        self.avl.cmd(' ')
        # --- set trim constraints ---
        if runTrim:
            CLreq = fc.get_CLreq()
            self.avl.cmd('A')
            self.avl.cmd('C')
            self.avl.cmd(CLreq)
            self.avl.cmd('D1') # elevator
            self.avl.cmd('PM')
            self.avl.cmd(fc.CmTrim)
        else:
            self.avl.cmd('A')
            self.avl.cmd('A')
            self.avl.cmd(alpha)
            self.avl.cmd('B')
            self.avl.cmd('B')
            self.avl.cmd(beta)
            self.avl.cmd('D1')
            self.avl.cmd('D1')
            self.avl.cmd(elevator)
        # --- run analysis ---
        self.avl.cmd('X')
        self.avl.cmd('ST')
        self.avl.cmd(' ')
        # --- system matrix ---
        #FIXME: dynamic output crashes AVL
#        self.avl.cmd(' ')
#        self.avl.cmd('mode')
#        self.avl.cmd('S')
        self.avl.cmd(' ')
        self.avl.cmd('QUIT')
    
    def _process_output(self):
        result = AVLresults(['e'])
        rawOutput = self.avl.get_output()
        coefficientsRaw=self._split_results(rawOutput,'Enter filename, or <return> for screen output   s>',' Neutral point  Xnp =',24)
        coefficientsRaw = coefficientsRaw.split()
        result = self._set_coefficients(coefficientsRaw,result)
        result = self._set_derivatives(coefficientsRaw,result)
        return result

    def _set_coefficients(self,raw,results):
        results.alpha= self._get_value(raw,'Alpha')
        results.beta =self._get_value(raw,'Beta')
        results.Mach =self._get_value(raw,'Mach')
        results.Sref =self._get_value(raw,'Sref')
        results.Cref =self._get_value(raw,'Cref')
        results.Bref =self._get_value(raw,'Bref')
        results.ARref=results.Bref**2/results.Sref
        results.coef.CX=self._get_value(raw,'CXtot')
        results.coef.CY=self._get_value(raw,'CYtot')
        results.coef.CZ=self._get_value(raw,'CZtot')
        results.coef.Cl=self._get_value(raw,'Cltot')
        results.coef.Cm=self._get_value(raw,'Cmtot')
        results.coef.Cn=self._get_value(raw,'Cntot')
        results.coef.CL=self._get_value(raw,'CLtot')
        results.coef.CD=self._get_value(raw,'CDtot')
        results.LD = results.coef.CL/results.coef.CD
        results.coef.CDind=self._get_value(raw,'CDind')
        results.k=results.coef.CDind/(results.coef.CL**2)
        results.e=1.0/(results.k*np.pi*results.ARref)       
        results.a=self._get_value(raw,'CLa')
        results.xNP=self._get_value(raw,'Xnp')
        results.CL0=results.coef.CL-results.a*np.pi/180*results.alpha
        results.elevator=self._get_value(raw,'elevator')
        return results

    def _set_derivatives(self,raw,results):
        results.derivs.CLa=self._get_value(raw,'CLa')
        results.derivs.CYa=self._get_value(raw,'CYa')
        results.derivs.Cla=self._get_value(raw,'Cla')
        results.derivs.Cma=self._get_value(raw,'Cma')
        results.derivs.Cna=self._get_value(raw,'Cna')
        
        results.derivs.CLb=self._get_value(raw,'CLb')
        results.derivs.CYb=self._get_value(raw,'CYb')         
        results.derivs.Clb=self._get_value(raw,'Clb')
        results.derivs.Cmb=self._get_value(raw,'Cmb')
        results.derivs.Cnb=self._get_value(raw,'Cnb')       
        
        results.derivs.CLp=self._get_value(raw,'CLp')
        results.derivs.CYp=self._get_value(raw,'CYp')
        results.derivs.Clp=self._get_value(raw,'Clp')
        results.derivs.Cmp=self._get_value(raw,'Cmp')
        results.derivs.Cnp=self._get_value(raw,'Cnp')
        
        results.derivs.CLp=self._get_value(raw,'CLq')
        results.derivs.CLq=self._get_value(raw,'CYq')
        results.derivs.CYq=self._get_value(raw,'Clq')
        results.derivs.Cmq=self._get_value(raw,'Cmq')
        results.derivs.Cnq=self._get_value(raw,'Cnq')

        results.derivs.CLr=self._get_value(raw,'CLr')
        results.derivs.CYr=self._get_value(raw,'CYr')
        results.derivs.Clr=self._get_value(raw,'Clr')
        results.derivs.Cmr=self._get_value(raw,'Cmr')
        results.derivs.Cnr=self._get_value(raw,'Cnr')

        results.derivs.CLde=self._get_value(raw,'CLd1')
        results.derivs.CYde=self._get_value(raw,'CYd1')
        results.derivs.Clde=self._get_value(raw,'Cld1')
        results.derivs.Cmde=self._get_value(raw,'Cmd1')
        results.derivs.Cnde=self._get_value(raw,'Cnd1')
        results.derivs.CDde=self._get_value(raw,'CDffd1')
        results.derivs.ede =self._get_value(raw,'ed1')
        return results

    def _split_results(self,source,fromStr1,toStr2,offset=0):
        a=self._findall(source,fromStr1)
        b=self._findall(source,toStr2)
        return source[a[1]:b[1]+offset]

    def _selTxt(self,source,fromStr1,toStr2=-1):
        a=self._findall(source,fromStr1)
        b=self._findall(source,toStr2)
        num=source[a[1]:b[0]]
        return float(num)

    def _findall(self,source,string):
        out = [(a.start(), a.end()) for a in list(re.finditer(string, source))]
        return out[0]
        
    def _selNum(self,source,fromStr,num):
        a=self._findall(source,fromStr)
        val=source[a[1]:a[1]+num]
        return float(val)

    def set_flight_conditions(self,flightConditions):
        pass
    
    def _get_value(self,raw,name,offset=2):
        idx = raw.index(name)
        val = raw[idx+2]
        return float(val)



class Aerodynamics(AVLsolver):
    def __init__(self,aircraft):
        super(Aerodynamics,self).__init__(aircraft)

class RunCases:
    def __init__(self):
        pass

def run_test1():
    import aircraft_FW as aircraft
    ac = aircraft.load('Baseline1')


    aero = Aerodynamics(ac)
    fc = FlightConditionsAVL(ac,0.7,1e4)
    results = aero.run_trim(fc)
    results.display()
    
    results2 = aero.run_single_point(fc,0.0,2.0)
    results2.display()


if __name__=="__main__":
    run_test1()