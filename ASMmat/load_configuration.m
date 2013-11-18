function missile = load_configuration(name)

if nargin<1
    load_configuration(1)
    return
end

missile.mass = 300;
missile.CG = [1.5, 0, 0];
missile.CD0 = 0.03;

missile.VLMpanels.wing = [6, 20];
missile.VLMpanels.tail = [6, 10];

missile.body.length = 2.5;
missile.body.diameter = 0.3;

missile.wing.secChords = [0.15, 0.1];
missile.wing.secIncidence = [2,1];
missile.wing.span = 3.0;
missile.wing.sweepAngle = 20;
missile.wing.dihedral = -5;
missile.wing.thickness = 0.12;
missile.wing.location = [0.6, 0.2];

missile.tail.secChords = [0.35, 0.15];
missile.tail.halfSpan = 0.4;
missile.tail.sweepAngle = 40;
missile.tail.thickness = 0.12;
missile.tail.locationX = 2.1;
missile.tail.centerOffset = 0.1;
missile.tail.elevatorChordRatio = 0.5;
missile.tail.numberOfTails = 4;
missile.tail.xAngle = 45;

end