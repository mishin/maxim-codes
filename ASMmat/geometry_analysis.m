function missile = geometry_analysis(missile)

if nargin<1
    missile = load_configuration(1);
    missile = geometry_analysis(missile);
    return
end

% all analysis for single segment wings only

%% general
missile.weight = missile.mass*9.81;

%% body
bodyFunction = fuselage_coordinate_generator(missile.body.type);
missile.body.coord = bodyFunction(missile.body.length,missile.body.diameter);
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
halfSpan = missile.fin.halfSpan;
missile.fin.areaPerSide = sum(missile.fin.secChords)*halfSpan/2;
missile.fin.area = missile.fin.areaPerSide*missile.fin.numberOfTails;

centerOffset = missile.fin.centerOffset;

secXm = centerOffset * tand(missile.fin.sweepAngle);
secYm = centerOffset;
secZm = 0.0;

secX = (halfSpan+centerOffset) * tand(missile.fin.sweepAngle);
secY = halfSpan+missile.fin.centerOffset;
secZ = 0.0;
rootSection = [0, 0, 0];
middleSection = [secXm, secYm,secZm];
tipSection  = [secX, secY, secZ];
missile.fin.rootSection = [];
missile.fin.middleSection = [];
missile.fin.tipSection = [];
centerOffsetRatio = missile.fin.secChords(1)/(missile.fin.secChords(1)+missile.fin.secChords(2));
missile.fin.middleChord = (missile.fin.secChords(2)+(missile.fin.secChords(1)-missile.fin.secChords(2))*(1-centerOffsetRatio));
for i=0:3
    angle = -missile.fin.xAngle - i*90;
    RotMatrix = [1, 0, 0; 0, cosd(angle), -sind(angle); 0, sind(angle), cosd(angle)];
    newRoot = rootSection*RotMatrix;
    newTip  = tipSection *RotMatrix;
    newMiddle = middleSection*RotMatrix;
    missile.fin.rootSection = [missile.fin.rootSection; newRoot];
    missile.fin.tipSection  = [missile.fin.tipSection; newTip];
    missile.fin.middleSection = [missile.fin.middleSection; newMiddle];
end

%% for mass analysis
missile.Mtot = missile.mass;

missile.fin.Xcg = missile.fin.locationX;
missile.fin.Nf = missile.fin.numberOfTails;
missile.fin.Area = missile.fin.areaPerSide;
missile.fin.AR = missile.fin.halfSpan^2/missile.fin.areaPerSide  *2;

missile.fuselage.l_tot = missile.body.length;
missile.fuselage.dmax = missile.body.diameter;
missile.fuselage.Xcg = 0.45*missile.body.length;

missile.Wtot_S = missile.mass/missile.wing.area;
missile.wing.tr = missile.wing.secChords(2) / missile.wing.secChords(1);
missile.wing.Cr = missile.wing.secChords(1);
missile.wing.Ct = missile.wing.secChords(2);
missile.wing.Area = missile.wing.area;
missile.wing.tc = missile.wing.thickness;
missile.wing.sweepLE = missile.wing.sweepAngle;
missile.wing.Xcg = missile.wing.location(1);
missile.wing.AR = missile.wing.span^2 / missile.wing.area;

end