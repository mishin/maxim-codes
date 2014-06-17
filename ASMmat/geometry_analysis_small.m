function config = geometry_analysis_small(config)

    function wing = get_PHIF(wing)
        step = 360/wing.numberOfTails;
        i = 1;
        for angle = 0:step:360-step
            wing.PHIF(i) = angle + wing.xAngle;
            i = i+1;
        end
    end

config.wing = get_PHIF(config.wing);
config.fin = get_PHIF(config.fin);


%% general
config.weight = config.mass*9.81;

%% body
bodyFunction = fuselage_coordinate_generator(config.body.type);
config.body.coord = bodyFunction(config.body.length,config.body.diameter);
config.body.sref = pi*config.body.diameter^2/4;
%% wing

config.wing.halfSpan = config.wing.span/2;
halfSpan = config.wing.halfSpan;
config.wing.area = sum(config.wing.secChords)*halfSpan;
secX = halfSpan*tand(config.wing.sweepAngle);
secY = halfSpan;
secZ = halfSpan*tand(config.wing.dihedral);
config.wing.rootSection = [0,0,0];
config.wing.tipSection = [secX, secY, secZ];

config.wing.MAC = mean(config.wing.secChords); % TODO: change MAC calculation
%% tail
halfSpan = config.fin.halfSpan;
config.fin.areaPerSide = sum(config.fin.secChords)*halfSpan/2;
config.fin.area = config.fin.areaPerSide*config.fin.numberOfTails;

centerOffset = config.fin.centerOffset;

secXm = centerOffset * tand(config.fin.sweepAngle);
secYm = centerOffset;
secZm = 0.0;

secX = (halfSpan+centerOffset) * tand(config.fin.sweepAngle);
secY = halfSpan+config.fin.centerOffset;
secZ = 0.0;
rootSection = [0, 0, 0];
middleSection = [secXm, secYm,secZm];
tipSection  = [secX, secY, secZ];
config.fin.rootSection = [];
config.fin.middleSection = [];
config.fin.tipSection = [];
centerOffsetRatio = config.fin.secChords(1)/(config.fin.secChords(1)+config.fin.secChords(2));
config.fin.middleChord = (config.fin.secChords(2)+(config.fin.secChords(1)-config.fin.secChords(2))*(1-centerOffsetRatio));
for i=0:3
    angle = -config.fin.xAngle - i*90;
    RotMatrix = [1, 0, 0; 0, cosd(angle), -sind(angle); 0, sind(angle), cosd(angle)];
    newRoot = rootSection*RotMatrix;
    newTip  = tipSection *RotMatrix;
    newMiddle = middleSection*RotMatrix;
    config.fin.rootSection = [config.fin.rootSection; newRoot];
    config.fin.tipSection  = [config.fin.tipSection; newTip];
    config.fin.middleSection = [config.fin.middleSection; newMiddle];
end

%% for mass analysis
config.Mtot = config.mass;

config.fin.Xcg = config.fin.locationX;
config.fin.Nf = config.fin.numberOfTails;
config.fin.Area = config.fin.areaPerSide;
config.fin.AR = config.fin.halfSpan^2/config.fin.areaPerSide  *2;

config.fuselage.l_tot = config.body.length;
config.fuselage.dmax = config.body.diameter;
config.fuselage.Xcg = 0.45*config.body.length;

config.Wtot_S = config.mass/config.wing.area;
config.wing.tr = config.wing.secChords(2) / config.wing.secChords(1);
config.wing.Cr = config.wing.secChords(1);
config.wing.Ct = config.wing.secChords(2);
config.wing.Area = config.wing.area;
config.wing.tc = config.wing.thickness;
config.wing.sweepLE = config.wing.sweepAngle;
config.wing.Xcg = config.wing.location(1);
config.wing.AR = config.wing.span^2 / config.wing.area;

end