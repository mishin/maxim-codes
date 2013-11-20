function config = load_configuration(name)

if nargin<1
    load_configuration(1)
    return
end

%mass.total = 300;
%mass.CG = [1.5, 0, 0];
%mass.inertia = [0,0,0];

config.mass = 300;
config.CG = [1.5 0 0];
config.CD0 = 0.03;

config.VLMpanels.wing = [6, 20];
config.VLMpanels.tail = [6, 10];

config.body.length = 2.5;
config.body.diameter = 0.3;

config.wing.secChords = [0.15, 0.1];
config.wing.secIncidence = [2,1];
config.wing.span = 3.0;
config.wing.sweepAngle = 20;
config.wing.dihedral = -5;
config.wing.thickness = 0.12;
config.wing.location = [0.6, 0.2];
config.wing.airfoil = 0;

config.tail.secChords = [0.35, 0.15];
config.tail.halfSpan = 0.4;
config.tail.sweepAngle = 40;
config.tail.thickness = 0.12;
config.tail.locationX = 2.1;
config.tail.centerOffset = 0.1;
config.tail.elevatorChordRatio = 0.5;
config.tail.numberOfTails = 4;
config.tail.xAngle = 45;

end