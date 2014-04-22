# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 21:29:18 2014

@author: Maxim
"""
import numpy as np
import constants
import quantities as qu
from weight_tools import MassList


class GeneralAviationMass:
    def __init__(self,aircraft,overrideColumn=0):
        self.aircraft = aircraft
        self.constants = constants.load('mass',overrideColumn)
        
        self.aircraftMass = aircraft.mass
        self.aircraftMass.airframe = MassList('Airframe')
        #atm                    =fc.ISAtmosphere(aircraft.designGoals.designAltitude)
        self._designQ          =aircraft.designGoals.fc.dynamicPressure
        self._designM          =aircraft.designGoals.fc.Mach
        self._Wdg              =float((aircraft.designGoals.designGrossMass*qu.kg).rescale('pound'))    
        self._Nz               =aircraft.designGoals.designLoadFactor*1.5
        self._q                =float((self._designQ*qu.Pa).rescale('pound_force/ft**2'))
        self._Nen              =float(len(aircraft.engine.engineList))

        self._wing()
        self._hStab()
        self._vStab()
        self._fuselage()
        self._landingGear()
        self._engine()
        self._fuel()
        self._misc()
#        self.emptyMass    =0.0
#        self.payloadMass  =0.0
#        self.totalMass    =0.0
        self.CG=np.zeros(3)
        self.get_inertia()
        self.aircraftMass.update_total()
    def _wing(self):
        #TODO: include equations here
        ac        =self.aircraft
        Sw        =float((ac.wing.area*qu.m**2).rescale('ft**2'))
        Wfw       =float((ac.designGoals.fuelMass*qu.kg).rescale('pound'))
        A         =ac.wing.aspectRatio                                      
        W_swept   =np.radians(ac.wing.sweep)
        taper     =ac.wing.taper
        t_c       =ac.wing.thicknessRatio
        Ww        =0.036*(Sw**0.758)*(Wfw**0.0035)*((A/((np.cos(W_swept))**2))**0.6)*(self._q**0.006)*((taper)**0.04)*(((100*(t_c ))/np.cos(W_swept))**(-0.3))*((self._Nz*self._Wdg)**0.49) 
        Mw        =float((Ww*qu.pound).rescale('kg'))
        CGR       =ac.wing.centroid
        CGR[0] = ac.wing.aapex[0]+0.4*ac.wing.MAC
        CGL       =np.array([1.,-1.,1.])*ac.wing.centroid
        inertia   =np.zeros(3)  #TODO - calculate inertia based on geometry
        if ac.wing.material=='composite':
            Mw *= 0.85
        self.aircraftMass.airframe.add_item('Left Wing',Mw/2.,CGL,inertia)
        self.aircraftMass.airframe.add_item('Right Wing',Mw/2.,CGR,inertia)

    def _hStab(self):
        ac        =self.aircraft
        A         =ac.hStab.aspectRatio
        Sht       =float((ac.hStab.area*qu.m**2).rescale('ft**2'))
        W_swept   =float(ac.hStab.sweep*qu.rad)
        taper_ht  =ac.hStab.taper
        t_c       =ac.hStab.thicknessRatio
        Wht       =0.016*((self._Nz*self._Wdg)**0.414)*(self._q**0.0168)*(Sht**0.896)*(((100*(t_c ))/np.cos(W_swept))**(-0.12))*((A/((np.cos(W_swept))**2))**0.043)*(taper_ht**(-0.02));
        Mht       =float((Wht*qu.pound).rescale('kg'))
        CGR       =ac.hStab.centroid
        CGL       =np.array([1.,-1.,1.])*ac.hStab.centroid
        inertia   =np.zeros(3)  #TODO - calculate inertia based on geometry
        if ac.hStab.material=='composite':
            Mht *= 0.83
        self.aircraftMass.airframe.add_item('Left Hstab',Mht/2.,CGL,inertia)
        self.aircraftMass.airframe.add_item('Right Hstab',Mht/2.,CGR,inertia)
    def _vStab(self):
        ac        =self.aircraft
        A         =ac.vStab.aspectRatio
        Svt       =float((ac.vStab.area*qu.m**2).rescale('ft**2'))
        W_swept   =float(ac.vStab.sweep*qu.rad)
        taper_vt  =ac.vStab.taper
        t_c       =ac.vStab.thicknessRatio
        Ht        =ac.hStab.aapex[2]-ac.vStab.aapex[2]
        Hv        =ac.vStab.aapex[2]+ac.vStab.span
        Ht = abs(Ht)
        Wvt       =0.073*(1.+(0.2*(Ht/Hv)))*((self._Nz*self._Wdg)**0.376)*(self._q**0.122)*(Svt**0.3873)*(((100*(t_c))/np.cos(W_swept))**(-0.49))*((A/((np.cos(W_swept))**2))**0.375)*(taper_vt**(0.039))
        Mvt       =float((Wvt*qu.pound).rescale('kg'))
        CG        =ac.vStab.centroid
        inertia   =np.zeros(3)  #TODO - calculate inertia based on geometry
        if ac.vStab.material=='composite':
            Mvt *= 0.83
        #Mvt *= 2.0 #TODO: find reliable mass calculation equation than this one
        Mvt *= 2.0
        self.aircraftMass.airframe.add_item('vStab',Mvt,CG,inertia)
    def _fuselage(self):
        #FIXME: fuselage mass prediction difference is due to wrong wetted area calculation
        ac        =self.aircraft
        Sf        =float((ac.fuselage.wettedArea*qu.m**2).rescale('ft**2'))
        xMW       =ac.wing.centroid[0]+ac.wing.MAC/4.
        xHT       =ac.hStab.centroid[0]+ac.hStab.MAC/4.
        Lt_si     =np.abs(xHT-xMW)
        Lt        =float((Lt_si*qu.m).rescale('ft'))
        L         =float((ac.fuselage.length*qu.m).rescale('ft'))
        D         =float((ac.fuselage.diameter*qu.m).rescale('ft'))
        Wpress    =0.
        Wfuselage =0.052*(Sf**1.086)*((self._Nz*self._Wdg)**0.177)*(Lt**(-0.051))*((L/D)**(-0.072))*(self._q**0.241)+Wpress
        Mfuselage =float((Wfuselage*qu.pound).rescale('kg'))
        if ac.fuselage.material=='composite':
            Wfuselage *= 0.9
        #CG        =ac.fuselage.centroid
        #TODO: check fuselage center of gravity prediction
        CG = self.constants.getValue('fuseCGratio') * self.aircraft.fuselage.length
        CG = np.array([CG,0.0,0.0])
        inertia   =np.zeros(3)
        self.aircraftMass.airframe.add_item('fuselage',Mfuselage,CG,inertia)
    def _landingGear(self):
        ac        =self.aircraft
        iNose     =np.argmin(ac.landingGear.groundContact_X)
        Nl        =ac.designGoals.designLandingLoadFactor*1.5        
        Wl        =self._Wdg        
        ln        =float((ac.landingGear.strutLength[iNose]*qu.m).rescale('ft'))
        WnoseGear =0.125*((Nl*Wl)**0.566)*((ln/12.)**0.845)
        MnoseGear =float((WnoseGear*qu.pound).rescale('kg'))
        CGnoseGear=np.array([ac.landingGear.groundContact_X[iNose],ac.landingGear.groundContact_Y[iNose],ac.landingGear.groundContact_Z[iNose]+ac.landingGear.tireDiameter[iNose]/2.])
        InoseGear =np.zeros(3)
        self.aircraftMass.airframe.add_item('nose gear',MnoseGear,CGnoseGear,InoseGear)
        nMG       =float(len(ac.landingGear.groundContact_X)-1.)
        for i in range(len(ac.landingGear.groundContact_X)):
            if i<>iNose:
                lm      =ac.landingGear.strutLength[i]
                Wmlg    =0.095*((Nl*Wl)**0.768)*((lm/12.)**0.409)
                Mmlg    =float((Wmlg*qu.pound).rescale('kg'))/nMG
                CG      =np.array([ac.landingGear.groundContact_X[i],ac.landingGear.groundContact_Y[i],ac.landingGear.groundContact_Z[i]+ac.landingGear.tireDiameter[i]/2.])
                inertia =np.zeros(3)
                name    ='landing gear '+str(i)
                self.aircraftMass.airframe.add_item(name,Mmlg,CG,inertia)
    def _engine(self):
        ac        =self.aircraft
        for i,engine in enumerate(ac.engine.engineList):
            Men  =engine.mass
            Wen  =float((Men*qu.kg).rescale('lb'))
            Weng =2.3*(Wen**0.922)
            Meng =float((Weng*qu.pound).rescale('kg'))
            CG   =np.array([ac.engine.CG_X[i],ac.engine.CG_Y[i],ac.engine.CG_Z[i]])
            I    =np.zeros(3)
            name ='engine '+str(i+1)
            self.aircraftMass.airframe.add_item(name,Meng,CG,I)
    def _fuel_(self,fuelFraction=1.0):
        # will be replaced, kept here for reference. Do not use!!!
        ac        =self.aircraft
        if ac.wing.fuelTank.mirrorFlag<>1:M=ac.wing.fuelTank.fuelMass*fuelFraction
        else:M=2.*ac.wing.fuelTank.fuelMass*fuelFraction
        CGR       =ac.wing.fuelTank.CG
        CGL       =CGR*np.array([1.,-1.,1])
        I         =np.zeros(3)
        self.aircraftMass.fuel.add_item('Fuel tank Right',M/2.,CGR,I)
        self.aircraftMass.fuel.add_item('Fuel tank Left' ,M/2.,CGL,I)
    def _fuel(self):
        sr = self.aircraft.wing.fuelTank.spanRatio
        cr = self.aircraft.wing.fuelTank.chordRatio
        CG = self.aircraft.wing.locateOnWing(sr,cr)
        self.aircraftMass.fuel.update_item_cg('Fuel tank Right',CG[0], CG[1],CG[2])
        self.aircraftMass.fuel.update_item_cg('Fuel tank Left', CG[0],-CG[1],CG[2])
    def _fuel_system(self):
        ac        =self.aircraft
        fuseCG    =(self.aircraftMass.airframe.get_item_by_name('fuselage')).coords
        I         =np.zeros(3)
        # Fuel System
        Vt        =float((ac.fuselage.diameter*qu.m).rescale('ft'))
        Vi        =0.0
        Nt        =2.0
        Wfs       =2.49*(Vt**0.726)*((1./(1.+(Vi/Vt)))**0.363)*(Nt**0.242)*(self._Nen**0.157)
        Mfs       =float((Wfs*qu.pound).rescale('kg'))
        name      ='fuel system'
        self.aircraftMass.airframe.add_item(name,Mfs,fuseCG,I)
    def _flight_controls(self):
        ac        =self.aircraft
        fuseCG    =(self.aircraftMass.airframe.get_item_by_name('fuselage')).coords
        I         =np.zeros(3)
        L         =float((ac.fuselage.length*qu.m).rescale('ft'))
        Bw        =float((ac.wing.span*qu.m).rescale('ft'))
        Wfc       =0.053*(L**0.1536)*(Bw**0.371)*((self._Nz*self._Wdg*(10.**(-4)))**0.8)
        Mfc       =float((Wfc*qu.pound).rescale('kg'))
        name      ='flight controls'
        self.aircraftMass.airframe.add_item(name,Mfc,fuseCG,I)
    def _hydraulics(self):
        ac        =self.aircraft
        fuseCG    =(self.aircraftMass.airframe.get_item_by_name('fuselage')).coords
        I         =np.zeros(3)
        Kh        =0.013
        W         =float((ac.fuselage.diameter*qu.m).rescale('ft'))
        M         =self._designM
        Whyd      =Kh*(W**0.8)*(M**0.5)
        Mhyd      =float((Whyd*qu.pound).rescale('kg'))
        name      ='hydraulics'
        self.aircraftMass.airframe.add_item(name,Mhyd,fuseCG,I)
    def _misc(self):
        self._fuel_system()
        self._flight_controls()
        self._hydraulics()
    def get_inertia(self):
        #Temporary statistical method until more detailed analysis becomes available
        #TODO: add detailed intertia analysis
        ac=self.aircraft 
        self.aircraftMass.airframe.update_totals()
        Mtot=self.aircraftMass.airframe.get_total_mass()
        R=[.25,.38,.39]
        Ixx=ac.wing.span**2*Mtot*R[0]**2/4.
        Iyy=ac.fuselage.length**2*Mtot*R[1]**2/4.
        Izz=(ac.wing.span/2+ac.fuselage.length/2)**2*Mtot*R[2]**2/4.
        self.inertia=np.array([Ixx,Iyy,Izz])
        self.aircraftMass.inertia = self.inertia