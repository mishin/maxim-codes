function config = load_configuration(name)

if nargin<1
    config = load_configuration(1);
    return
end

config.type = 1; % 1 - high AR, 2 - small AR
config.mass = 93;
%config.CG = [0.92 0 0];
config.Vmax = 640/3.6; %m/sec
config.CD0 = 0.03;
config.liftToWeigthRatio = 0.5;

config.VLMpanels.wing = [6, 20];
config.VLMpanels.tail = [8, 8];

config.body.length = 1.76;
config.body.diameter = 0.18;
config.body.type = 'cylinder';

config.wing.secChords = [0.131, 0.1025];
config.wing.secIncidence = [0,0];
config.wing.span = 1.68;
config.wing.sweepAngle = 13; % leading edge sweep
config.wing.dihedral = 0;
config.wing.thickness = 0.15;
config.wing.location = [0.6838, 0.05]; % XZ
%config.wing.airfoil = 0;
config.wing.skinThickness = 0.0005;
config.wing.airfoil = 'NASALRN1015.txt';

config.fin.secChords = [0.08, 0.08];
config.fin.halfSpan = 0.1;
config.fin.sweepAngle = 0;
config.fin.thickness = 0.12;
config.fin.locationX = 1.6868;
config.fin.centerOffset = config.body.diameter/2;
config.fin.elevatorChordRatio = 0;
config.fin.numberOfTails = 4;
config.fin.xAngle = 45;
config.fin.airfoil = 'NACA0012.txt';

config = geometry_analysis(config);
config.massData = weight_analysis(config);
config.CG = [config.massData.Xcg 0 0];

end