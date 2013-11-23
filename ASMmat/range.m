function [R,aero] = range(gb, hStart, hEnd, Vstart)

if nargin==0
    clc
    gb = load_configuration(1);
    tic
    [R,~] = range(gb,8500, 1000, 640/3.6);
    toc
    return
end

% assumption: constant velocity at lift = 0.6 Weight

nSeg = 10;
hStep = (hEnd-hStart)/nSeg;

V = ones(nSeg,1) * Vstart;
dR = zeros(nSeg,1);
h = hStart:hStep:hEnd;
dt= zeros(nSeg,1);
LD = zeros(nSeg,1);

V(1) = Vstart;
for i = 1:nSeg;
    alt = h(i);
    fc = get_flight_conditions(alt,V(i));
    aero = missile_trim_AVL(gb, fc);
    dR(i) = -hStep*aero.LD;
    dt(i) = dR(i)/V(i);
%    alpha(i) = aero.alpha;
%    de(i) = aero.elevator;
%    LD(i) = aero.LD;
%    SM(i) = aero.SM;
end
R =sum(dR)/1e3;
t = sum(dt)/60;

end