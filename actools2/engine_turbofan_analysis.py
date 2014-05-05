# -*- coding: utf-8 -*-
"""
Created on Fri May 02 15:37:24 2014

@author: Maxim
"""
from numpy import exp,pi,sqrt, zeros, arange

#function [d02, thrustD, sfcD, sfcH]= EngineModeling(threq, Mcrz, AltH, thr, hh, Mf)
def engine_modeling(threq, Mcrz, AltH, thr, hh, Mf):

    """
    threq : SLS Design thrust 
    Mcrz : Design Cruise Mach number
    AltH : Design Cruise Altitude [m]
    thr : thrust required for Flight Condition 
    hh  : Flight Altitude
    Mf  : Flight Mach Number
    
    d02 : Fan diameter
    thrustD : thrust at Design Cruise Condition
    sfcD : sfc at Design Cruise Condition
    sfcF : sfc at Flight Condition
    """
     
# engine design data
    bprdes = 0.8
    pifan = 2.345
    etafan = 0.84
    pilpc = 3.0
    etalpc = 0.84
    pihpc = 6.50
    etahpc = 0.84
    pib = 0.94
    titdes = 1600
    etab = 0.95
    etahpt = 0.88
    etalpt = 0.90 
    etanz = 0.98
    m20  = 0.5
    dspinner = 0.20
    pimx = 0.97
    
    
# Sea Level Static Condition Design
    
    tt2, pt2,_1,_2 = airproperty(0)

# Core flow LPC    
    tt23, pt23, worklpc=comp(pilpc, etalpc, tt2, pt2)
# Fan flow      
    tt13, pt13, workfan=comp(pifan, etafan, tt2, pt2)
# Core flow HPC   
    tt3, pt3, workhpc=comp(pihpc, etahpc, tt23, pt23)
    
    
# Burner Calculation
    tt4 = titdes
    pt4, ff, phi, a40=burner(tt4, tt3, pt3, pib, etab)
    
    
    
# HPT Calculation    
    work = workhpc 
    tt45, pt45=turbine(etahpt, tt4, pt4, phi, work)
    
# LPT Calculation    
    work = worklpc + bprdes*workfan 
    tt5, pt5=turbine(etalpt, tt45, pt45, phi, work)
    

# Mixer Calculation
    ttmx, ptmx, phimx= Mixer(bprdes, pt5, pt13, tt5, tt13, pimx, phi) 

# Sonic Nozzle Calculation   of mixed flow

    p6, vjet, area6= nozzle(etanz, ttmx, ptmx, phimx)
    
# Thrust Calculation   
    t0, p0,_1,_2 = airproperty(0)
    thrust = ( vjet + (p6-p0)*area6 )*(1+bprdes)
    
# sfc Calculation       
    sfc = (ff*3600*9.8)/(thrust) 
    
# Engine Diameter   
    d02, mfr, mcore= EngineSize(threq, thrust, bprdes, dspinner, m20, tt2,pt2) 
    
# Engine Lip Area Design
    a1   = 0.25*pi*( d02**2 - dspinner**2) 
    a2   = a1
    a4   = a40*mcore     # NGV Area
    
     
# Cycle analsys specific to M0, AltH, thrH      
     
     
     
     
# Design Point at 65000 ft (19500m)

    pid    = 0.98  
     
    t0,p0,rho0,viscos=airproperty(AltH)   
    gamma, Cp0 = cal_constant(t0, 0)
    u0        = Mcrz*sqrt(gamma*287.058*t0) 
 
    fact1   = ( 1 + .5*(gamma-1)*Mcrz**2 ) 
    fact2   = fact1**( gamma/(gamma-1))
    
    tt0 = t0*fact1 
    pt0 = p0*fact2
    
     
# Intake : Adiabatic with 2# loss
   
    tt2   = tt0  
    pt2  = pid*pt0 
        
# Core flow LPC   
    tt23, pt23, worklpc=comp(pilpc, etalpc, tt2, pt2)
        
# Fan flow      
    tt13, pt13, workfan=comp(pifan, etafan, tt2, pt2)
    
    
# Core flow HPC   
    tt3, pt3, workhpc=comp(pihpc, etahpc, tt23, pt23)   
    
    
# Burner Calculation
    tt4 = titdes
    pt4, ff, phi, a40=burner(tt4, tt3, pt3, pib, etab)
    
    
# HPT Calculation    
    work = workhpc 
    tt45, pt45=turbine(etahpt, tt4, pt4, phi, work)
    
# LPT Calculation    
    work = worklpc + bprdes*workfan 
    tt5, pt5=turbine(etalpt, tt45, pt45, phi, work)
    

# Mixer Calculation
    ttmx, ptmx, phimx= Mixer(bprdes, pt5, pt13, tt5, tt13, pimx, phi) 

# Sonic Nozzle Calculation   of mixed flow

    p6, vjet, area6= nozzle(etanz, ttmx, ptmx, phimx)
    

    
# Thrust Calculation   

    thrust = ( vjet - u0 + (p6-p0)*area6 )*(1+bprdes)
    
# Mass flow rate Calculation    
    mfr, mcore= MassFlow(bprdes, a2, m20, tt2,pt2)     

    Dadd, Dvis= InletDrag(AltH, Mcrz , a1, mfr) 
     
    NDrag= NozzleDrag(AltH, Mcrz , d02)  
    
    
    thrust = thrust*mcore - (Dadd + Dvis + NDrag) 
     
    thrustD = thrust 
    
    Anze    = area6*mfr 

# sfc Calculation       
    sfc = (mcore*ff*3600*9.8)/(thrust) 
 
    
    sfcD = sfc
         
# Operation Line     
    gamma2, Cp2 = cal_constant(tt2, 0)
    Mnorm2= normass(gamma2, m20) 
    gamma4, Cp4 = cal_constant(tt4, phi)
    Mnorm4= normass(gamma4, 1.0) 
    
    CT = sqrt(tt2/tt4)
    CP = sqrt(Cp2/Cp4)
    CA = (a4*(1+bprdes)/a2) 
    CM = 1/(1+ff) 

    Cfac = CM*CT*CA*pib*Mnorm4 
    
    TR  = tt4/tt2 
   
    M21 = 0.3 
    M22 = 0.55 
    
    for i in range(1,21):
        mm  = 0.5*(M22 + M21)
        Mnorm2 = normass(1.4,mm) 
        pic   = Mnorm2/Cfac 

        pifan1  = pic*pifan/(pilpc*pihpc) 
        pilpc1  = pic/(pihpc)   
        pihpc1 = pic/(pilpc) 

        pid    = 0.98  

        t0,p0,rho0,viscos=airproperty(hh) 
        
        gamma, Cp0 = cal_constant(t0, 0)
        
        u0        = Mf*sqrt(gamma*287.058*t0) 

        fact1   = ( 1 + .5*(gamma-1)*Mf**2 ) 
        fact2   = fact1**( gamma/(gamma-1))    

        tt0 = t0*fact1 
        pt0 = p0*fact2

        # Intake : Adiabatic with 2# loss

        tt2   = tt0
        pt2  = pid*pt0 

        # Core flow LPC   
        tt23, pt23, worklpc=comp(pilpc1, etalpc, tt2, pt2)

        # Fan flow      
        tt13, pt13, workfan=comp(pifan1, etafan, tt2, pt2)

        # Core flow HPC   
        tt3, pt3, workhpc=comp(pihpc1, etahpc, tt23, pt23)   

        # Burner Calculation
        tt4 = TR*tt2 
        if tt4 > titdes:
            tt4 = titdes

        pt4, ff, phi, a40=burner(tt4, tt3, pt3, pib, etab)

        # HPT Calculation    
        work = workhpc 
        tt45, pt45=turbine(etahpt, tt4, pt4, phi, work)

        # LPT Calculation    
        work = worklpc + bprdes*workfan 
        tt5, pt5=turbine(etalpt, tt45, pt45, phi, work)


        # Mixer Calculation
        ttmx, ptmx, phimx= Mixer(bprdes, pt5, pt13, tt5, tt13, pimx, phi) 

        # Sonic Nozzle Calculation   of mixed flow

        p6, vjet, area6= nozzle(etanz, ttmx, ptmx, phimx)
        

        # Thrust Calculation   

        thrust = ( vjet - u0 + (p6-p0)*area6 )*(1+bprdes)    
           

        # Mass flow rate Calculation    
        mfr, mcore= MassFlow(bprdes, a2, mm, tt2,pt2)     

        Dadd, Dvis= InletDrag(hh, Mf , a1, mfr) 

        NDrag= NozzleDrag(hh, Mf , d02)  


        thrust = thrust*mcore - (Dadd + Dvis + NDrag) 

        # sfc Calculation       
        sfc = (mcore*ff*3600*9.8)/(thrust)
        
        if thrust>thr:
            M22 = mm 
        else:
            M21 = mm

    sfcH = sfc 
    Anze2 = area6*mfr
    tt44  = tt4
    thrustt = thrust
    mach = mm

    return d02, thrustD, sfcD, sfcH




def NozzleDrag(AltH, M0, d02):
# Function calculates nozzle base drag

    gamma = 1.4 
    Cp = gamma*(287.058)/(gamma-1.0) 

    t0,p0,rho0,viscos=airproperty(AltH) 
   
    u0        = M0*sqrt(gamma*287.058*t0) 
    
    CD       = 0.02 
    
    NDrag   = CD*( 0.5*rho0*u0**2*(0.25*pi*d02**2.0) ) 
    return NDrag


#function [Dadd, Dvis]= InletDrag(AltH, M0, a1, mfr)
def InletDrag(AltH, M0, a1, mfr):
# Function calculates additive inlet drag

    gamma = 1.4 
    Cp = gamma*(287.058)/(gamma-1.0) 
    fact1   = ( 1.0 + .5*(gamma-1.0)*M0**2 ) 
    fact2   = fact1**( gamma/(gamma-1.0))    
    
    t0,p0,rho0,viscos=airproperty(AltH) 
    
    tt0 = t0*fact1 
    pt0 = p0*fact2
    
    Mnorm = normass(gamma, M0) 
    
    area2   = sqrt(Cp*tt0)/(Mnorm*pt0)   
    a0       = area2*mfr  
    
    M1  = 0.0 
    M2  = 1.0 
    
    for i in range(1,21):
    
        mm = 0.5*(M1 + M2) 
        Mnorm= normass(gamma, mm) 

        area2   = sqrt(Cp*tt0)/(Mnorm*pt0)   
        a01      = area2*mfr

        dela = a1 - a01 
        
        if dela>0.0:
            M2    = mm 
        else:
            M1    = mm 

    
    fact1   = ( 1 + .5*(gamma-1)*mm**2 ) 
    fact2   = fact1**( gamma/(gamma-1))    

    p1       = pt0/fact2 
    t1        = tt0/fact1 
    u1        = mm*sqrt(gamma*287.058*t1) 
    rho1     = p1/(287.058*t1) 
    
    Dadd = p1*a1*(1 + gamma*mm**2) - p0*a0*( (a1/a0) + gamma*M0**2 ) 
    
    r1     = sqrt(a1/pi) 
    ld     = 8*r1 
   
    Rel  = rho1*u1*ld/viscos 
   
    cfl   = 0.074*( Rel**(-0.2) )  
    
    Dvis = cfl*( 0.5*rho1*u1**2 )*(2*pi*ld*r1)
    
    return Dadd, Dvis




#function [ttmx, ptmx, phimx]= Mixer(bpr, pt5, pt15, tt5, tt15, pimx, phi)
def Mixer(bpr, pt5, pt15, tt5, tt15, pimx, phi):
# Function calculates engine size

    gamma1, Cp1 = cal_constant(tt15,0)
    gamma2, Cp2 = cal_constant(tt5,phi)

    ptmx = pimx*(pt5 + bpr*pt15)/(1 + bpr) 
    
    ttmx = (tt5 + bpr*tt15)/(1 + bpr) 
    phimx = (phi + bpr)/(1+bpr) 
    return ttmx, ptmx, phimx

#function [mfr, mcore]= MassFlow(bpr, a2, m20, tt2,pt2)
def MassFlow(bpr, a2, m20, tt2,pt2):
# Function calculates mass flow rates

    gamma = 1.4 
    Cp = gamma*(287.058)/(gamma-1) 

    Mnorm= normass(gamma, m20) 
    
    area2 = sqrt(Cp*tt2)/(Mnorm*pt2) 
    
    mfr = a2/area2      
    mcore = mfr/(1+bpr)  
    return mfr, mcore




#function [d02, mfr, mcore]= EngineSize(threq, thmfr, bpr, dspinner, m20, tt2,pt2)
def EngineSize(threq, thmfr, bpr, dspinner, m20, tt2,pt2):
# Function calculates engine size

    gamma = 1.4 
    Cp = gamma*(287.058)/(gamma-1) 

    Mnorm= normass(gamma, m20) 
    
    area2 = sqrt(Cp*tt2)/(Mnorm*pt2) 
    
    mcore = threq/thmfr  
    mfr = (1+bpr)*mcore   
    a02   = area2*mfr + 0.25*pi*dspinner**2 
    d02    = 2*sqrt(a02/pi) 
    return d02, mfr, mcore


#function [tt3, pt3, work]=comp(pic, etac, tt2, pt2)
# Function calculates compressor performance
def comp(pic, etac, tt2, pt2):

    temp = tt2
    phi = 0
    [gamma11, Cp11] = cal_constant(temp,phi)
    gamma22 = gamma11
    Cp22 = Cp11
    
    for i in range(1,21): 
        gamma = 0.5*(gamma11 + gamma22)
        Cpav = 0.5*(Cp11 + Cp22)
        gamma1 = (gamma-1)/(gamma)
        tt3 = tt2*( (pic**gamma1 -1)/etac + 1 )
        [gamma22, Cp22] = cal_constant(tt3,phi)
    
    work = Cpav*(tt3-tt2)
    pt3 = pt2*pic
    return tt3, pt3, work


#function [pt4, ff, phi, a4]=burner(tt4, tt3, pt3, pib, etab)
def burner(tt4, tt3, pt3, pib, etab):
# Function calculates burner performance

    hfl = 43.124*(10**6)
    pt4 = pt3*pib
    phi = 0
    [gamma, Cp1] = cal_constant(tt3,phi)
    phi = 0.5
    
    for i in range(1,21):
        gamma, Cp2 = cal_constant(tt4,phi)
        ff = 0.5*(Cp1 +Cp2)*(tt4-tt3)/(etab*hfl)
        phi = ff/0.0676
   
    Mnorm= normass(gamma, 1.0)
    
    a4 = sqrt(Cp2*tt4)/(Mnorm*pt4)
    return pt4, ff, phi, a4



#function [tt5, pt5]=turbine(etat, tt4, pt4, phi, work)
# Function calculates turbine performance
def turbine(etat, tt4, pt4, phi, work):
    temp = tt4
    [gamma11, Cp11] = cal_constant(temp,phi)
    gamma22 = gamma11
    Cp22 = Cp11
    
    for i in range(1,21):
        Cpav = 0.5*(Cp11 + Cp22)
        gamma = 0.5*(gamma11 + gamma22)
        tt5 = tt4 - work/Cpav 
        taut = tt5/tt4
        gamma10 = gamma/(gamma-1)
        pit = ( (taut -1)/etat + 1 )**(gamma10)
        gamma22, Cp22 = cal_constant(tt5,phi)

    pt5 = pt4*pit
    return tt5, pt5

#function [p6, vjet, area6]= nozzle(etanz, tt5, pt5, phi)
# Function calculates Sonic Nozzle performance
def nozzle(etanz, tt5, pt5, phi):
    mach   = 1.0
    [gamma11, Cp11] = cal_constant(tt5,phi)
    gamma22 = gamma11
    Cp22 = Cp11
    
    for i in range(1,21):
        gamma = 0.5*(gamma11 + gamma22)
        Cpav = 0.5*(Cp11 + Cp22)
        fact1 = ( 1 + .5*(gamma -1)*mach**2 ) 
        t6is = tt5/fact1 
        t6 = tt5 - etanz*(tt5-t6is) 
        gamma22, Cp22 = cal_constant(t6,phi)
    
    vjet = sqrt(287.058*gamma22*t6) 
    
    tt6    = tt5
    fact2 = t6/t6is 
    
    fact3 = ( 1- (1-fact2)/etanz )**(gamma22/(gamma22-1))  
    fact4 = ( 1 + .5*(gamma22 -1)*mach**2 )**(gamma22/(gamma22-1))  
    pin = fact2**(gamma22/(gamma22-1))  
    
    pt6 = pt5/pin 
    p6 = pt6/fact4 
   
    Mnorm= normass(gamma22, mach)
    
    area6 = sqrt(Cp22*tt6)/(Mnorm*pt6) 
    return p6, vjet, area6


def normass(gamma, mach):
# Function calculates normalized mass flow rate

   gamm1 = gamma-1 
   gamp1 = gamma+1 
   Mnorm = mach*( gamma/sqrt(gamm1) )*( (1 + .5*(gamm1)*mach**2 )**(-.5*gamp1/gamm1) ) 
   return Mnorm



def cal_constant(temp,equivratio):

# Function calculates lpc performance
    phi     = equivratio
    Tgamma  = 3055.5555556
    RR      = (287.058)
    
    gamma   = 1 + (0.4)/( 1 + 0.4*((Tgamma/temp)**2 * exp(Tgamma/temp) /(exp(Tgamma/temp) -1)**2)) 
    gamma   = gamma - (0.004184*phi + 0.000096*phi**2)  
    Cpgas   = RR*(gamma)/(gamma-1.)
    return gamma,Cpgas

    

def airproperty(altitude):
# Function produces standard air property data for a given altitude in [m]
# temp in [K]
# press in [Pa]
# density in [kg/m**3]

    if altitude <11000.:
        tempt = 288.15 - 0.00650*altitude 
        press = 101325*( (tempt/288.15)**(5.2569) )
    else:
        tempt = 216.65 
        press = 22632*exp(1.73-0.0001573*altitude)


    density = press/(286.7*tempt)
    c1 =  (1.458e-6)
    ss = 110.5
    
    viscos = c1*(tempt**1.5)/(tempt + ss)
    
    return tempt,press,density,viscos


def run_test1():
    import matplotlib.pyplot as plt
    Treq = 7.9e3
    McrD = 0.9
    altD = 0.0
    TreqF = 5.0e3
    Mach = arange(0.1,1.1,0.1)
    alt = arange(0.0,15e3,2e3)
    T = zeros([len(alt),len(Mach)])
    plt.figure(1)
    plt.grid(True)
    plt.hold(True)
    plt.xlabel('Mach number')
    plt.ylabel('SFC')
    legend = list()
    for i, h in enumerate(alt):
        for j, M in enumerate(Mach):
            d, Td, sfcD, sfcF = engine_modeling(Treq,McrD,altD,TreqF,h,M)
            T[i,j] = Td
        legend.append('h=%.0fm'%h)
        plt.plot(Mach,T[i])
    plt.legend(legend,'upper left')
    plt.show()
    
    #print EngineModeling(Treq,McrD,AltD,TreqF,h,Mf)
#d02, thrustD, sfcD, sfcH
if __name__=="__main__":
    run_test1()