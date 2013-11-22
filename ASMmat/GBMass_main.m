% temporary main program for Mass Analysis 
% programed by Kwon-Su Jeon

clc
clear all

Config.Wtot_S = 747.4741;               % wing loading (kgf/m2)   
Config.Mtot = 93;                       % total Mass (kg)  


% wing configuration
Config.wing.tr = 1;                     % wing taper ratio (non-dim.)
Config.wing.Cr = 0.0889;                % wing root chord (m)   
Config.wing.Ct = Config.wing.Cr*Config.wing.tr;              % tip chord length (m)
Config.wing.Area = 0.1244;              % root chord thickness ratio
Config.wing.tc = 0.15;                  % root chord thickness ratio
Config.wing.sweepLE = 0;                % wing LE sweep angle (deg.)
Config.wing.span = 1.3995;              % span length (m)   
Config.wing.skinThickness = 0.0005;          % Skin thickness (m)  ! 두께가 작아질 수록 중량 증가?
%Config.wing.AR = 19.2;              % wing AR (non-dim.)
Config.wing.AR = 15.74241;              % wing AR (non-dim.)
Config.wing.Xcg = 0.9;                  % wing x-location (m)


% fuselage configuration
Config.fuselage.l_tot = 1.7784;         % fuselage structural length (m)
Config.fuselage.dmax = 0.1778;          % fuselage maximum diameter (m)
Config.fuselage.Xcg = 0.45*Config.fuselage.l_tot;     % fuselage x-location (m)
Config.body.length = Config.fuselage.l_tot;

% fin configuration
Config.fin.Area = 0.00762;              % exposed area per 1 fin (unit : m^2)
Config.fin.AR = 2.5;                    % fin AR for one fin
Config.fin.Nf = 4;                      % no. of fin
Config.fin.Xcg = 1.6;                   % fin x-location (m)

%payload
%Config.payload.Xcg = 0.5*Config.fuselage.l_tot;               % wing x-location (m)

Vmax = 640;                             % max. velocity (km/h)

%input for mass analysis

MassIn.payloadFusRatio = 0.5;
MassIn.fuselageCGratio = 0.45;
MassIn.Vmax = Vmax;                     % max. velocity (km/h)
MassIn.Xcg_seeker = 0.04425;            % x-location of seeker (m)   
MassIn.Xcg_imu = 0.3135;                % x-location of IMU (m)
MassIn.Xcg_gps = 0.3135;                % x-location of GPS (m)
MassIn.Xcg_seekele = 0.5385;            % x-location of seeker electronics (m)
MassIn.Xcg_misscom = 1.35624;           % x-location of mission computer (m)
MassIn.Xcg_batt = 1.47882;              % x-location of thermal battery (m)
MassIn.Xcg_cas = 1.6899;                % x-location of control actuator system (m)

% call Mass analysis part
%geom = load_configuration(1);
[MassOut] = GBMass(Config, MassIn);
