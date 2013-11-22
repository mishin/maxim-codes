% Mass Estimation Module
% programed by Kwon-Su Jeon
%
function [MassOut]=GBMass(Config, MassIn) 

format short

%% unit conversion factor

unit = 1;
U_m2ft = convlength(unit,'m','ft');    % unit conversion factor (m -> ft)
U_m2in = convlength(unit,'m','in');    % unit conversion factor (m -> in)
U_lb2kg = convmass(unit,'lbm','kg');    % unit conversion factor (lbs -> kg) 
U_kg2lb = convmass(unit,'kg','lbm');    % unit conversion factor (kg -> lbs) 
Ukph2fps = 0.911344;       % unit conversion factor (km/h  -> ft/sec)      

% parameters
g=32.2;                    % gravity acceleration in Britich Unit (ft/sec^2)
n = 1;                     % loading factor (non-dim.) 
faert = 0.75;              % aeroelastic tailoring factor  
fcomp = 0.75;              % composite material factor 
fult = 3.5;                % ultimate load factor
Kvs = 1.19;                % factor for variable sweep angle wing (from Raymer's text book p.592 ) 
%fm_f = 4;                  % material factor for fin, weight per area : 4 lbs/ft2 

Mpod = 0;                           % Engine pod mass (kg)   
Sflap = 0;                          % flap area (m2)   

%% component weight estimation
% wing Weight : large type (high AR case)
tr = Config.wing.tr;                % wing taper ratio (non-dim.)
Cr = Config.wing.Cr;                % wing root chord (m)   
Ct = Cr*tr;                         % tip chord length (m)
tc = Config.wing.tc;                % root chord thickness ratio
Sw = Config.wing.Area;              % wing area (m2)

Wtot_S = Config.Wtot_S;             % wing loading (kgf/m2)   
ARw = Config.wing.AR;               % wing AR (non-dim.)
Mtot = Config.Mtot;                 % total Mass (kg)   
SweepW = Config.wing.sweepLE;       % wing LE sweep angle (deg.)

bw = Config.wing.span;              % span length (m)   
t_r=Config.wing.skinThickness;          %   Skin thickness (m)

l_tot = Config.fuselage.l_tot;      % fuselage structural length (m)
dmax = Config.fuselage.dmax;        % fuselage maximum diameter (m)

t_r =  convlength(t_r,'m','ft');

Aexp_f = Config.fin.Area;           % exposed area per 1 fin (unit : m^2)
AR_f = Config.fin.AR;               % fin AR for one fin
Nf = Config.fin.Nf;                 % no. of fin

Vmax = MassIn.Vmax;                 % max. velocity (km/h)
Xcg_w = Config.wing.Xcg;            % wing cg x-location (m)
Xcg_fus = Config.fuselage.Xcg;      % fuselage cg x-location (m)
Xcg_fin = Config.fin.Xcg;           % fin cg x-location (m)

% subsystem X-location
Xcg_seeker = MassIn.Xcg_seeker; 
Xcg_imu = MassIn.Xcg_imu;
Xcg_gps = MassIn.Xcg_gps;
Xcg_seekele = MassIn.Xcg_seekele;
Xcg_misscom = MassIn.Xcg_misscom;
Xcg_batt = MassIn.Xcg_batt;
Xcg_cas = MassIn.Xcg_cas;

% payload X-location
Config.payload.Xcg = Config.body.length*MassIn.payloadFusRatio;
Xcg_pl = Config.payload.Xcg;

%Sw = Mtot*g/(Wtot_S);            % wing area (m2)   
sbw=bw/2;                  % semi-span length (m)
sbw_ft = convlength(sbw,'m','ft');
bw_ft = convlength(bw,'m','ft');
Cr_ft = convlength(Cr,'m','ft');
Ct_ft = convlength(Ct,'m','ft');


Mtot_lb = convmass(Mtot,'kg','lbm');
%% wing mass estimation using FLOPS method

syms x y z
Cy=Cr_ft-2*y*Cr_ft*(1-tr)/bw_ft;

% Load distribution over the wing
q_lw=2*Mtot_lb*g*n*(sbw_ft^2-x^2)^0.5/sbw_ft^2/pi;
%q_lw=2*Mtot_lb*n*(sbw_ft^2-x^2)^0.5/sbw_ft^2/pi;
% Bending moment is assumed only the lift components
MBend=int(q_lw*x,0,y);

% The flange areas
Sflange=MBend/(t_r*Cy*cos(SweepW*pi/180));      % flange area (ft^2)

% The volume required
Vol=vpa(int(Sflange,0,1),5);

% The total lift distribution over the wing
Lift=vpa(int(q_lw,0,1),5);

% The bending material factor calculation
if (ARw>5)
    Ka=ARw-5;
else
    Ka=0;
end
d = ARw^0.25*(1+(0.5*faert-0.16)*(sin(SweepW*pi/180))^2+0.3*Ka*(1-0.5*faert)*sin(SweepW*pi/180));
Bt = 2*Vol/(Lift*d);
Bte = 0;                    % Due to the absence of the engine on the wing

%%%%% Calculate the wing weight based on FLOPS
Ke = 1-(Bte/Bt)*(Mpod/Mtot);
K = 8.8*10^(-6)*(1+sqrt(6.25/bw_ft))*Bt;

Wtot_lb = convmass(Mtot,'kg','lbm');                        % Total weight (kg -> lbs)

Sw_ft = Sw*U_m2ft^2;                                        % Wing area (m2 -> ft2)

Wwing1 = K*fult*bw_ft*(1-0.4*fcomp)*(1-0.1*faert);          % wing bending material weight (lbs)
Wwing2 = 0.68*(1-0.17*fcomp)*Sflap^0.34*Wtot_lb^0.6;        % wing shear material and flaps weight (lbs)
Wwing3 = 0.035*(1-0.3*fcomp)*Sw_ft^1.5 ;                    % wing contraol surfaces and non-structural weight (lbs)    

% The wing weight estimation
Wwing = (Wtot_lb*Ke*Wwing1+Wwing2+Wwing3)/(1+Wwing1)*Kvs;    % wing weight (lbs)

Mw = double(Wwing);
Mw = convmass(Mw,'lbm','kg');

% wing weight : small type (low AR case)

% fin weight
Af = U_m2ft^2*Aexp_f;       % exposed wing area, unit conversion (m2 to ft2)
AR_AS = AR_f;       % fin AR
E_AS = U_m2ft^2*Aexp_f;  % exposed wing area, unit conversion (m2 to ft2)

Mfin_p1_lb = 6.77483 *(E_AS)^1.02 *(AR_AS)^0.56;       % fin weight per each panel (lbs) 

%Mfin_p1_lb = fm_f*Af;                                 % fin weight per each panel (lbs)  

Mfin_lb = Mfin_p1_lb*Nf;                               % total fin weight (lbs) 

Mfin_p1 = convmass(Mfin_p1_lb,'lbm','kg');             % fin weight per each panel (kg)
Mfin = convmass(Mfin_lb,'lbm','kg');                % total fin weight (kg)


% fuselage weight
Lbs = l_tot* U_m2in;                       % fuselage structural length (m -> in)
Dbs = dmax* U_m2in;                        % fuselage maximum diameter (m -> in)

Mbs = 0.0604*(Lbs^0.64) * (Dbs^1.77);      % fuselage structural weight (lbs)
Mbs_lb = Mbs;
Mbs = Mbs_lb*U_lb2kg;                      % total fuselage weight (kg)

%% system weight
% nose part 
M_seeker = 0;                           % seeker weight (kg)

% center body 1 part 
M_gps = 0.125;                          % GPS weight (kg)
M_imu = 0.24;                           % IMU weight (kg)
M_seekele = 0;                          % seeker electronics weight (kg)

% center body 2 part
M_misscom = 0.906;                      % misison computer weight (kg)
M_batt = 2.718;                         % thermal battery (kg)

% aft body part

Vmax = Vmax*Ukph2fps;                   % max velocity (ft/sec)
Sf = 0.061;                             % fin exposed area (unit : m^2)
Sf = Sf*U_m2ft^2;
Nfcs = 1;                               % number of flight control system

M_cas = 0.00002*Vmax^2*Sf*Nfcs;         % control actuation system weight (lbs)
M_cas = convmass(M_cas,'lbm','kg');     % control actuation system weight (kg)

Msub = M_seeker + M_gps + M_imu + M_seekele + M_misscom + M_batt + M_cas;

%% airframe Mass 
Mairframe = Mw + Mfin + Mbs + Msub;     % airframe mass (kg)

Mpl = Mtot - Mairframe;                 % payload mass (kg)


%% AC cg estimation

co1 = Xcg_w*Mw + Xcg_fus*Mbs + Xcg_fin*Mfin;
co2 = Xcg_seeker*M_seeker + Xcg_imu*M_imu + Xcg_gps*M_gps + Xcg_seekele*M_seekele;
co3 = Xcg_misscom*M_misscom + Xcg_batt*M_batt + Xcg_cas*M_cas + Xcg_pl*Mpl;

Xcg = (co1+ co2 +co3)/Mtot;


%% Mass analysis results


MassOut.Mw = Mw;                        % wing mass (kg)
MassOut.Mfin_p1 = Mfin_p1;              % fin mass per 1 fin (kg)
MassOut.Mfin = Mfin;                    % total fin mass (kg)
MassOut.Mbs = Mbs;                      % body structural mass (kg) : fuselage
MassOut.Msub = Msub;                    % subsystem mass (kg)
MassOut.Mairframe = Mairframe;          % airframe mass (kg)
MassOut.Mpl = Mpl;                      % payload mass (kg)
MassOut.Xcg = Xcg;                      % total cg x-location (m)


end




