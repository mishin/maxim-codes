function config = load_configuration(name)

if nargin<1
    load_configuration(1)
    return
end

config.mass = 93;
config.CG = [0.92 0 0];
config.CD0 = 0.03;
config.liftToWeigthRatio = 0.5;

config.VLMpanels.wing = [6, 20];
config.VLMpanels.tail = [8, 8];

config.body.length = 2.5;
config.body.diameter = 0.3;

config.wing.secChords = [0.0889, 0.0889];
config.wing.secIncidence = [1,1];
config.wing.span = 1.3995;
config.wing.sweepAngle = 0;
config.wing.dihedral = -5;
config.wing.thickness = 0.12;
config.wing.location = [0.9, 0.0];
config.wing.airfoil = 0;

config.tail.secChords = [0.0762, 0.0762];
config.tail.halfSpan = 0.1;
config.tail.sweepAngle = 0;
config.tail.thickness = 0.12;
config.tail.locationX = 1.6;
config.tail.centerOffset = 0.01;
config.tail.elevatorChordRatio = 0;
config.tail.numberOfTails = 4;
config.tail.xAngle = 45;

config = geometry_analysis(config);

end