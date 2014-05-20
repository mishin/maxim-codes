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
            CD0 = aircraft.get_drag()
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
        pth.set_file_prefix('flyingwing')
        path = pth.get_tmp_file('avl')
        cg = self.ac.get_cg()
        CD0 = 0.0 #self.ac.get_drag()
        fid = open(path,'wt')
        # general input
        fid.write('%s\n'%self.ac.name)
        fid.write('0.0\n')
        fid.write('0  0  0.0\n')
        fid.write('%.4f  %.4f  %.4f\n'%(self.ac.wing.area,self.ac.wing.MAC,self.ac.wing.span))
        fid.write('%.4f  %.4f  %.4f\n'%(cg[0],cg[1],cg[2]))
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
            fid.write('%.6f  %.6f  %.6f  %.6f  %.6f\n'%(apex[0],apex[1],apex[2],chord,angle))
            fid.write('AFIL\n%s\n'%pathAf)
            fid.write('CONTROL\n')
            fid.write('elevator  1.0  %.4f  0  0  0  +1\n'%(self.ac.wing.elevon.location[i]))
        fid.close()
    
    def run_trim(self,fc):
        self.create_input_file()
        self.avl = AVL()
        self.run_avl_trim(fc)
        result = self._process_output()
        result.CD0 = fc.CD0
        result.SM = float((result.xNP - fc.cg[0])/fc.Lref)
        self.avl.terminate()
        return result
        #self.avl.terminate()

    def run_avl_trim(self,fc):
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
        CLreq = fc.get_CLreq()
        self.avl.cmd('A')
        self.avl.cmd('C')
        self.avl.cmd(CLreq)
        self.avl.cmd('D1') # elevator
        self.avl.cmd('PM')
        self.avl.cmd(fc.CmTrim)
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
        coefficientsRaw=self._split_results(rawOutput,'Enter filename, or <return> for screen output   s>',' Neutral point  Xnp =',12)      
        result = self._set_coefficients(coefficientsRaw,result)
        result = self._set_derivatives(coefficientsRaw,result)
        return result

    def _set_coefficients(self,raw,results):
        results.alpha=self._selTxt(raw,'Alpha =','pb/2V')
        results.beta =self._selTxt(raw,'Beta  =','qc/2V')
        results.Mach =self._selTxt(raw,'Mach  =','rb/2V')
        results.Sref =self._selTxt(raw,'Sref =','Cref =')
        results.Cref =self._selTxt(raw,'Cref =','Bref =')
        results.Bref =self._selTxt(raw,'Bref =','    \n  Xref')
        results.ARref=results.Bref**2/results.Sref
        results.coef.CX=self._selTxt(raw,'CXtot =','Cltot =')
        results.coef.CY=self._selTxt(raw,'CYtot =','Cmtot =')
        results.coef.CZ=self._selTxt(raw,'CZtot =','Cntot =')
        results.coef.Cl=self._selTxt(raw,'Cltot =','Cl\'tot =')
        results.coef.Cm=self._selTxt(raw,'Cmtot =','\n  CZtot =')
        results.coef.Cn=self._selTxt(raw,'Cntot =','Cn\'tot =')
        results.coef.CL=self._selTxt(raw,'CLtot =','\n  CDtot')
        results.coef.CD=self._selTxt(raw,'CDtot =','\n  CDvis')
        results.coef.CDind=self._selTxt(raw,'CDff  =','\| Trefftz\n  CYff ')
        results.k=results.coef.CDind/(results.coef.CL**2)         
        results.e=1.0/(results.k*np.pi*results.ARref)       
        results.a=self._selTxt(raw,'CLa =','CLb =')
        results.xNP=self._selNum(raw,'Neutral point  Xnp =   ',10)
        results.CL0=results.coef.CL-results.a*np.pi/180*results.alpha
        results.elevator=self._selNum(raw,'elevator        =',10)
        return results

    def _set_derivatives(self,raw,results):
        results.derivs.CLa=self._selTxt(raw,'CLa =','CLb =')            
        results.derivs.CYa=self._selTxt(raw,'CYa =','CYb =')
        results.derivs.Cla=self._selTxt(raw,'Cla =','Clb =')
        results.derivs.Cma=self._selTxt(raw,'Cma =','Cmb =')
        results.derivs.Cna=self._selTxt(raw,'Cna =','Cnb =')
        
        results.derivs.CLb=self._selNum(raw,'CLb =',11)
        results.derivs.CYb=self._selNum(raw,'CYb =',11)            
        results.derivs.Clb=self._selNum(raw,'Clb =',11)
        results.derivs.Cmb=self._selNum(raw,'Cmb =',11)
        results.derivs.Cnb=self._selNum(raw,'Cnb =',11)            
        
        results.derivs.CLp=self._selTxt(raw,'CLp =','CLq =')      
        results.derivs.CYp=self._selTxt(raw,'CYp =','CYq =')      
        results.derivs.Clp=self._selTxt(raw,'Clp =','Clq =')      
        results.derivs.Cmp=self._selTxt(raw,'Cmp =','Cmq =')      
        results.derivs.Cnp=self._selTxt(raw,'Cnp =','Cnq =')                  
        
        results.derivs.CLp=self._selTxt(raw,'CLq =','CLr =')  
        results.derivs.CLq=self._selTxt(raw,'CYq =','CYr =')  
        results.derivs.CYq=self._selTxt(raw,'Clq =','Clr =')  
        results.derivs.Cmq=self._selTxt(raw,'Cmq =','Cmr =')  
        results.derivs.Cnq=self._selTxt(raw,'Cnq =','Cnr =')  

        results.derivs.CLr=self._selNum(raw,'CLr =',11)
        results.derivs.CYr=self._selNum(raw,'CYr =',11)
        results.derivs.Clr=self._selNum(raw,'Clr =',11)
        results.derivs.Cmr=self._selNum(raw,'Cmr =',11)
        results.derivs.Cnr=self._selNum(raw,'Cnr =',11)

        results.derivs.CLde=self._selNum(raw,'CLd1 =',11)
        results.derivs.CYde=self._selNum(raw,'CYd1 =',11)
        results.derivs.Clde=self._selNum(raw,'Cld1 =',11)
        results.derivs.Cmde=self._selNum(raw,'Cmd1 =',11)
        results.derivs.Cnde=self._selNum(raw,'Cnd1 =',11)
        results.derivs.CDfe=self._selNum(raw,'CDffd1 =',11)
        results.derivs.ede =self._selNum(raw,'ed1 =',11)
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
    def run_single_point(self):
        pass
    def run_alpha_sweep(self):
        pass
    def run_beta_sweep(self):
        pass


class Aerodynamics(AVLsolver):
    def __init__(self,aircraft):
        super(Aerodynamics,self).__init__(aircraft)

class RunCases:
    def __init__(self):
        pass

def run_test1():
    from misc_tools import Timer
    import aircraft_FW as aircraft
    timer = Timer()
    ac = aircraft.load('X47B')
    ac.display()
    timer.lap('load')
    #ac.display()
    
    aero = Aerodynamics(ac)
    fc = FlightConditionsAVL(ac,200,1e4)
    aero.run_trim(fc).display()
    timer.stop('aero')
if __name__=="__main__":
    run_test1()