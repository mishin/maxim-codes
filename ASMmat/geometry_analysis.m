function missile = geometry_analysis(missile)

if nargin<1
    missile = load_configuration(1);
    missile = geometry_analysis(missile);
    return
end

% all analysis for single segment wings only

%% general
missile.weight = missile.mass*9.81;

%% wing

missile.wing.halfSpan = missile.wing.span/2;
halfSpan = missile.wing.halfSpan;
missile.wing.area = sum(missile.wing.secChords)*halfSpan;
secX = halfSpan*tand(missile.wing.sweepAngle);
secY = halfSpan;
secZ = halfSpan*tand(missile.wing.dihedral);
missile.wing.rootSection = [0,0,0];
missile.wing.tipSection = [secX, secY, secZ];

missile.wing.MAC = mean(missile.wing.secChords); % TODO: change MAC calculation
%% tail
halfSpan = missile.tail.halfSpan;
missile.tail.areaPerSide = sum(missile.tail.secChords)*halfSpan/2;
missile.tail.area = missile.tail.areaPerSide*missile.tail.numberOfTails;

centerOffset = missile.tail.centerOffset;
secX = halfSpan * tand(missile.tail.sweepAngle);
secY = halfSpan + centerOffset;
secZ = 0.0;
rootSection = [0, centerOffset, 0];
tipSection  = [secX, secY, secZ];
missile.tail.rootSection = [];
missile.tail.tipSection = [];

for i=0:3
    angle = -missile.tail.xAngle - i*90;
    RotMatrix = [1, 0, 0; 0, cosd(angle), -sind(angle); 0, sind(angle), cosd(angle)];
    newRoot = rootSection*RotMatrix;
    newTip  = tipSection *RotMatrix;
    missile.tail.rootSection = [missile.tail.rootSection; newRoot];
    missile.tail.tipSection  = [missile.tail.tipSection; newTip];
end

end