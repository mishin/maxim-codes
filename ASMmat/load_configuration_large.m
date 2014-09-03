function config = load_configuration_large()

config.type = 1; % 1- large, 2- small
config.mass = 93; %kg
config.Vmax = 640/3.6; %m/sec

config.CD0 = 0.03; % assumed value

config.VLMpanels.wing = [6, 20];
config.VLMpanels.tail = [8, 8];

config.body.length = 1.76;
config.body.diameter = 0.18;
config.body.type = 'cylinder';

config.wing.secChords = [0.131, 0.131];
config.wing.secIncidence = [0,0];
config.wing.span = 1.68;
config.wing.sweepAngle = 13; % leading edge sweep
config.wing.dihedral = -3;
config.wing.location = [0.6838, 0.05]; % XZ
config.wing.skinThickness = 0.0005; % for weight estimation
config.wing.thickness = 0.15;
config.wing.airfoil = 'NASALRN1015.txt';

config.fin.secChords = [0.08, 0.08];
config.fin.halfSpan = 0.1;
config.fin.sweepAngle = 0;
config.fin.locationX = 1.6868;
config.fin.centerOffset = 0.9*config.body.diameter/2; % first section offset
config.fin.elevatorChordRatio = 0; %[0;1] 0 - all moving
config.fin.numberOfTails = 4;
config.fin.xAngle = 45;
config.fin.thickness = 0.12;
config.fin.airfoil = 'NACA0012.txt';

config = geometry_analysis_large(config);
config.massData = weight_analysis_large(config);
config.CG = [config.massData.Xcg 0 0];
config.drag = drag_analysis(config);

end