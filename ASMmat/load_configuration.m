function config = load_configuration(name)

if nargin<1
    config = load_configuration(1);
    return
end

config.mass = 93;
%config.CG = [0.92 0 0];
config.Vmax = 640/3.6; %m/sec
config.CD0 = 0.03;
config.liftToWeigthRatio = 0.5;

config.VLMpanels.wing = [6, 20];
config.VLMpanels.tail = [8, 8];

config.body.length = 1.7784;
config.body.diameter = 0.1778;

config.wing.secChords = [0.0889, 0.0889];
config.wing.secIncidence = [1,1];
config.wing.span = 1.3995;
config.wing.sweepAngle = 0;
config.wing.dihedral = 0;
config.wing.thickness = 0.15;
config.wing.location = [0.9, 0.2];
config.wing.airfoil = 0;
config.wing.skinThickness = 0.0005;

config.fin.secChords = [0.0762, 0.0762];
config.fin.halfSpan = 0.1;
config.fin.sweepAngle = 0;
config.fin.thickness = 0.12;
config.fin.locationX = 1.6;
config.fin.centerOffsetRatio = 0.1;
config.fin.elevatorChordRatio = 0;
config.fin.numberOfTails = 4;
config.fin.xAngle = 45;

config = geometry_analysis(config);
config.massData = weight_analysis(config);
config.CG = [config.massData.Xcg 0 0];

end